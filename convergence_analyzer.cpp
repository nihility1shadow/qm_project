#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kTiny = 1.0e-15;

struct Options {
  std::vector<std::string> data_files;
  std::string manifest;
  std::vector<std::string> metadata_files;
  std::map<std::string, std::string> parameters;
  std::map<std::string, double> predict_parameters;
  std::string value_columns = "4:";
  std::vector<std::string> series_names;
  std::string out_prefix = "convergence";
  int time_column = 0;
  int segment_count = 5;
  int subblocks = 4;
  int min_tail_segments = 2;
  int max_plot_points = 1600;
  double segment_width = 0.0;
  double absolute_tolerance = 1.0e-4;
  double relative_tolerance = 0.01;
  double pass_probability = 0.95;
  double ridge_penalty = 1.0;
  bool self_test = false;
};

struct DataSet {
  std::string path;
  std::vector<double> time;
  std::vector<std::vector<double>> series;
};

struct Metric {
  double estimate = 0.0;
  double standard_error = 0.0;
  double probability = 0.0;
};

struct SeriesSegment {
  double mean = 0.0;
  double mean_se = 0.0;
  double noise_sigma = 0.0;
  double rho1 = 0.0;
  Metric amplitude;
  Metric drift;
  Metric jump;
  double score = 0.0;
};

struct SegmentResult {
  double start = 0.0;
  double end = 0.0;
  std::size_t begin_index = 0;
  std::size_t end_index = 0;
  std::size_t point_count = 0;
  std::vector<SeriesSegment> series;
  double score = 0.0;
  double tail_score = 0.0;
  int worst_series = -1;
  std::string limiting_metric;
  std::string classification;
};

struct AnalysisResult {
  std::string case_id;
  std::map<std::string, std::string> parameters;
  std::vector<std::string> data_files;
  std::vector<std::string> series_names;
  std::vector<double> time;
  std::vector<std::vector<double>> mean_series;
  std::vector<double> noise_rms_time;
  std::vector<double> tolerances;
  std::vector<SegmentResult> segments;
  double total_index = 0.0;
  double best_eligible_tail_score = 0.0;
  double confirmed_start = std::numeric_limits<double>::quiet_NaN();
  double possible_start = std::numeric_limits<double>::quiet_NaN();
  double noise_growth = std::numeric_limits<double>::quiet_NaN();
  std::string status;
  std::string warning;
};

struct ManifestCase {
  std::string id;
  std::vector<std::string> files;
  std::map<std::string, std::string> parameters;
};

struct TendencyRow {
  std::string parameter;
  double mean = 0.0;
  double scale = 0.0;
  double coefficient = 0.0;
  std::string direction;
};

struct TendencyModel {
  bool fitted = false;
  bool reliable = false;
  std::vector<TendencyRow> rows;
  std::vector<double> coefficients;
  double intercept = 0.0;
  double rmse = 0.0;
  double predicted_index = std::numeric_limits<double>::quiet_NaN();
  double predicted_low = std::numeric_limits<double>::quiet_NaN();
  double predicted_high = std::numeric_limits<double>::quiet_NaN();
  std::string note;
};

std::string trim(const std::string &text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string &text, char delimiter) {
  std::vector<std::string> parts;
  std::stringstream stream(text);
  std::string part;
  while (std::getline(stream, part, delimiter)) parts.push_back(trim(part));
  return parts;
}

double parse_double(const std::string &text, const std::string &label) {
  char *end = nullptr;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(value)) {
    throw std::runtime_error("invalid " + label + ": " + text);
  }
  return value;
}

int parse_int(const std::string &text, const std::string &label) {
  const double value = parse_double(text, label);
  if (std::floor(value) != value || value < 0.0 ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("invalid " + label + ": " + text);
  }
  return static_cast<int>(value);
}

std::pair<std::string, std::string> parse_key_value(const std::string &text,
                                                    const std::string &label) {
  const std::size_t equals = text.find('=');
  if (equals == std::string::npos || equals == 0 || equals + 1 >= text.size()) {
    throw std::runtime_error(label + " must use NAME=VALUE: " + text);
  }
  return {trim(text.substr(0, equals)), trim(text.substr(equals + 1))};
}

std::string parent_path(const std::string &path) {
  const std::size_t pos = path.find_last_of("/\\");
  return pos == std::string::npos ? "" : path.substr(0, pos);
}

bool is_absolute_path(const std::string &path) {
  return (!path.empty() && (path[0] == '/' || path[0] == '\\')) ||
         (path.size() > 1 && path[1] == ':');
}

std::string join_path(const std::string &base, const std::string &path) {
  if (base.empty() || is_absolute_path(path)) return path;
  const char separator = base.find('\\') != std::string::npos ? '\\' : '/';
  return base + separator + path;
}

std::string html_escape(const std::string &text) {
  std::string out;
  out.reserve(text.size());
  for (char c : text) {
    if (c == '&') out += "&amp;";
    else if (c == '<') out += "&lt;";
    else if (c == '>') out += "&gt;";
    else if (c == '"') out += "&quot;";
    else out += c;
  }
  return out;
}

std::string json_escape(const std::string &text) {
  std::string out;
  out.reserve(text.size() + 8);
  for (char c : text) {
    switch (c) {
      case '\\': out += "\\\\"; break;
      case '"': out += "\\\""; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default: out += c;
    }
  }
  return out;
}

std::string csv_escape(const std::string &text) {
  if (text.find_first_of(",\"\r\n") == std::string::npos) return text;
  std::string out = "\"";
  for (char c : text) out += c == '"' ? "\"\"" : std::string(1, c);
  return out + "\"";
}

std::vector<std::string> parse_csv_row(const std::string &line) {
  std::vector<std::string> cells;
  std::string cell;
  bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (quoted) {
      if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') {
        cell += '"';
        ++i;
      } else if (c == '"') {
        quoted = false;
      } else {
        cell += c;
      }
    } else if (c == '"') {
      quoted = true;
    } else if (c == ',') {
      cells.push_back(trim(cell));
      cell.clear();
    } else {
      cell += c;
    }
  }
  cells.push_back(trim(cell));
  return cells;
}

std::vector<double> parse_numeric_row(const std::string &line) {
  std::vector<double> values;
  const char *cursor = line.c_str();
  while (*cursor != '\0') {
    while (*cursor != '\0' &&
           (std::isspace(static_cast<unsigned char>(*cursor)) || *cursor == ',' ||
            *cursor == ';')) {
      ++cursor;
    }
    if (*cursor == '\0' || *cursor == '#') break;
    char *end = nullptr;
    errno = 0;
    const double value = std::strtod(cursor, &end);
    if (end == cursor || errno == ERANGE || !std::isfinite(value)) {
      throw std::runtime_error("non-numeric token near: " +
                               std::string(cursor, std::min<std::size_t>(30, std::strlen(cursor))));
    }
    values.push_back(value);
    cursor = end;
  }
  return values;
}

std::vector<int> resolve_columns(const std::string &spec, int column_count,
                                 int time_column) {
  std::vector<int> columns;
  for (const std::string &token : split(spec, ',')) {
    if (token.empty()) continue;
    const std::size_t colon = token.find(':');
    if (colon == std::string::npos) {
      columns.push_back(parse_int(token, "value column"));
    } else {
      const std::string left = trim(token.substr(0, colon));
      const std::string right = trim(token.substr(colon + 1));
      const int first = left.empty() ? 0 : parse_int(left, "column range start");
      const int last = right.empty() ? column_count - 1
                                     : parse_int(right, "column range end");
      if (last < first) throw std::runtime_error("column range is reversed: " + token);
      for (int column = first; column <= last; ++column) columns.push_back(column);
    }
  }
  std::sort(columns.begin(), columns.end());
  columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
  columns.erase(std::remove(columns.begin(), columns.end(), time_column), columns.end());
  if (columns.empty()) throw std::runtime_error("no value columns selected");
  for (int column : columns) {
    if (column < 0 || column >= column_count) {
      throw std::runtime_error("value column " + std::to_string(column) +
                               " is outside data width " +
                               std::to_string(column_count));
    }
  }
  return columns;
}

DataSet read_data(const std::string &path, int time_column,
                  const std::string &column_spec) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open data file: " + path);
  DataSet result;
  result.path = path;
  std::vector<int> columns;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::string cleaned = trim(line);
    if (cleaned.empty() || cleaned[0] == '#') continue;
    std::vector<double> row;
    try {
      row = parse_numeric_row(cleaned);
    } catch (const std::exception &error) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) + ": " +
                               error.what());
    }
    if (row.empty()) continue;
    if (time_column < 0 || time_column >= static_cast<int>(row.size())) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) +
                               ": time column is outside the row");
    }
    if (columns.empty()) {
      columns = resolve_columns(column_spec, static_cast<int>(row.size()), time_column);
      result.series.assign(columns.size(), {});
    }
    const int required = *std::max_element(columns.begin(), columns.end());
    if (required >= static_cast<int>(row.size())) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) +
                               ": inconsistent column count");
    }
    const double time = row[time_column];
    if (!result.time.empty() && !(time > result.time.back())) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) +
                               ": time must be strictly increasing");
    }
    result.time.push_back(time);
    for (std::size_t j = 0; j < columns.size(); ++j) {
      result.series[j].push_back(row[columns[j]]);
    }
  }
  if (result.time.size() < 12) {
    throw std::runtime_error(path + ": at least 12 numeric rows are required");
  }
  return result;
}

void validate_repeats(const std::vector<DataSet> &sets) {
  if (sets.empty()) throw std::runtime_error("at least one data file is required");
  const DataSet &base = sets.front();
  for (std::size_t r = 1; r < sets.size(); ++r) {
    if (sets[r].time.size() != base.time.size() ||
        sets[r].series.size() != base.series.size()) {
      throw std::runtime_error("repeat files have different shapes: " + sets[r].path);
    }
    for (std::size_t i = 0; i < base.time.size(); ++i) {
      const double scale = 1.0 + std::fabs(base.time[i]);
      if (std::fabs(sets[r].time[i] - base.time[i]) > 1.0e-9 * scale) {
        throw std::runtime_error("repeat files have different time grids: " + sets[r].path);
      }
    }
  }
}

bool identical_observations(const DataSet &left, const DataSet &right) {
  if (left.time != right.time || left.series.size() != right.series.size()) {
    return false;
  }
  for (std::size_t j = 0; j < left.series.size(); ++j) {
    if (left.series[j] != right.series[j]) return false;
  }
  return true;
}

double normal_cdf(double value) {
  return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

double sample_quantile(std::vector<double> values, double probability) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double position = probability * static_cast<double>(values.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(position));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - static_cast<double>(lo);
  return values[lo] * (1.0 - fraction) + values[hi] * fraction;
}

double robust_range(const std::vector<double> &values) {
  const std::size_t maximum = 10000;
  std::vector<double> sample;
  if (values.size() <= maximum) {
    sample = values;
  } else {
    sample.reserve(maximum);
    for (std::size_t i = 0; i < maximum; ++i) {
      const std::size_t index = i * (values.size() - 1) / (maximum - 1);
      sample.push_back(values[index]);
    }
  }
  return std::max(0.0, sample_quantile(sample, 0.95) - sample_quantile(sample, 0.05));
}

double mean_range(const std::vector<double> &values, std::size_t begin,
                  std::size_t end) {
  if (end <= begin) return 0.0;
  return std::accumulate(values.begin() + static_cast<std::ptrdiff_t>(begin),
                         values.begin() + static_cast<std::ptrdiff_t>(end), 0.0) /
         static_cast<double>(end - begin);
}

std::pair<double, double> linear_fit(const std::vector<double> &time,
                                     const std::vector<double> &value,
                                     std::size_t begin, std::size_t end) {
  const double mean_t = mean_range(time, begin, end);
  const double mean_y = mean_range(value, begin, end);
  double cross = 0.0;
  double square = 0.0;
  for (std::size_t i = begin; i < end; ++i) {
    const double centered = time[i] - mean_t;
    cross += centered * (value[i] - mean_y);
    square += centered * centered;
  }
  const double slope = square > kTiny ? cross / square : 0.0;
  return {mean_y - slope * mean_t, slope};
}

double residual_sigma(const std::vector<double> &time,
                      const std::vector<double> &value, std::size_t begin,
                      std::size_t end, double intercept, double slope) {
  std::vector<double> residuals;
  residuals.reserve(end - begin);
  for (std::size_t i = begin; i < end; ++i) {
    residuals.push_back(value[i] - intercept - slope * time[i]);
  }
  const double median = sample_quantile(residuals, 0.5);
  for (double &residual : residuals) residual = std::fabs(residual - median);
  const double mad = sample_quantile(residuals, 0.5);
  return std::max(kTiny, mad / 0.6744897501960817);
}

double lag_one_correlation(const std::vector<std::vector<double>> &residuals,
                           std::size_t begin, std::size_t end) {
  if (end <= begin + 2) return 0.0;
  double numerator = 0.0;
  double denominator_previous = 0.0;
  double denominator_current = 0.0;
  for (const auto &residual : residuals) {
    for (std::size_t i = begin + 1; i < end; ++i) {
      numerator += residual[i - 1] * residual[i];
      denominator_previous += residual[i - 1] * residual[i - 1];
      denominator_current += residual[i] * residual[i];
    }
  }
  const double denominator = std::sqrt(denominator_previous * denominator_current);
  if (denominator <= kTiny) return 0.0;
  return std::max(-0.5, std::min(0.95, numerator / denominator));
}

double effective_count(double count, double rho) {
  const double result = count * (1.0 - rho) / std::max(kTiny, 1.0 + rho);
  return std::max(2.0, std::min(count, result));
}

double small_sample_inflation(int degrees_of_freedom) {
  static const double t95[] = {
      0.0, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
      2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131,
      2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069,
      2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};
  if (degrees_of_freedom <= 0) return 3.0;
  if (degrees_of_freedom <= 30) return t95[degrees_of_freedom] / 1.959963984540054;
  return 1.0;
}

Metric make_metric(double estimate, double standard_error, double tolerance,
                   double inflation) {
  Metric metric;
  metric.estimate = std::max(0.0, estimate);
  metric.standard_error = std::max(kTiny, standard_error * inflation);
  metric.probability = normal_cdf((tolerance - metric.estimate) /
                                  metric.standard_error);
  metric.probability = std::max(0.0, std::min(1.0, metric.probability));
  return metric;
}

std::string classify(double score) {
  if (score >= 95.0) return "strong";
  if (score >= 80.0) return "usable";
  if (score >= 50.0) return "uncertain";
  return "not_converged";
}

std::vector<std::pair<std::size_t, std::size_t>> make_segments(
    const std::vector<double> &time, int segment_count, double segment_width) {
  std::vector<std::pair<std::size_t, std::size_t>> segments;
  const double first = time.front();
  const double last = time.back();
  int count = segment_count;
  if (segment_width > 0.0) {
    count = std::max(1, static_cast<int>(std::ceil((last - first) / segment_width)));
  }
  for (int segment = 0; segment < count; ++segment) {
    const double start = first + (last - first) * segment / count;
    const double end = first + (last - first) * (segment + 1) / count;
    const auto begin_it = std::lower_bound(time.begin(), time.end(), start);
    const auto end_it = segment + 1 == count
                            ? time.end()
                            : std::lower_bound(time.begin(), time.end(), end);
    const std::size_t begin = static_cast<std::size_t>(begin_it - time.begin());
    const std::size_t finish = static_cast<std::size_t>(end_it - time.begin());
    if (finish > begin + 3) segments.push_back({begin, finish});
  }
  if (segments.empty()) throw std::runtime_error("segments contain too few points");
  return segments;
}

AnalysisResult analyze_case(const std::string &case_id,
                            const std::vector<std::string> &files,
                            const std::map<std::string, std::string> &parameters,
                            const Options &options) {
  std::vector<DataSet> sets;
  sets.reserve(files.size());
  for (const std::string &file : files) {
    sets.push_back(read_data(file, options.time_column, options.value_columns));
  }
  validate_repeats(sets);
  std::vector<DataSet> unique_sets;
  std::vector<std::string> duplicate_paths;
  unique_sets.reserve(sets.size());
  for (DataSet &set : sets) {
    const auto duplicate = std::find_if(
        unique_sets.begin(), unique_sets.end(),
        [&set](const DataSet &candidate) {
          return identical_observations(set, candidate);
        });
    if (duplicate == unique_sets.end()) {
      unique_sets.push_back(std::move(set));
    } else {
      duplicate_paths.push_back(set.path);
    }
  }
  sets = std::move(unique_sets);
  const std::size_t repeat_count = sets.size();
  const std::size_t point_count = sets[0].time.size();
  const std::size_t series_count = sets[0].series.size();

  AnalysisResult result;
  result.case_id = case_id;
  result.parameters = parameters;
  for (const DataSet &set : sets) result.data_files.push_back(set.path);
  result.time = sets[0].time;
  result.mean_series.assign(series_count, std::vector<double>(point_count, 0.0));
  result.noise_rms_time.assign(point_count, 0.0);
  result.series_names = options.series_names;
  if (result.series_names.empty()) {
    for (std::size_t j = 0; j < series_count; ++j) {
      result.series_names.push_back("series_" + std::to_string(j));
    }
  }
  if (result.series_names.size() != series_count) {
    throw std::runtime_error("series name count does not match selected columns");
  }

  for (std::size_t j = 0; j < series_count; ++j) {
    for (std::size_t i = 0; i < point_count; ++i) {
      for (const DataSet &set : sets) result.mean_series[j][i] += set.series[j][i];
      result.mean_series[j][i] /= static_cast<double>(repeat_count);
    }
  }
  if (repeat_count >= 2) {
    for (std::size_t i = 0; i < point_count; ++i) {
      double sum = 0.0;
      for (std::size_t j = 0; j < series_count; ++j) {
        double variance = 0.0;
        for (const DataSet &set : sets) {
          const double delta = set.series[j][i] - result.mean_series[j][i];
          variance += delta * delta;
        }
        variance /= static_cast<double>(repeat_count - 1);
        sum += variance;
      }
      result.noise_rms_time[i] = std::sqrt(sum / static_cast<double>(series_count));
    }
  }

  result.tolerances.resize(series_count);
  for (std::size_t j = 0; j < series_count; ++j) {
    result.tolerances[j] = options.absolute_tolerance +
                           options.relative_tolerance * robust_range(result.mean_series[j]);
  }

  const auto ranges = make_segments(result.time, options.segment_count,
                                    options.segment_width);
  result.segments.resize(ranges.size());
  for (std::size_t k = 0; k < ranges.size(); ++k) {
    SegmentResult &segment = result.segments[k];
    segment.begin_index = ranges[k].first;
    segment.end_index = ranges[k].second;
    segment.point_count = segment.end_index - segment.begin_index;
    segment.start = result.time[segment.begin_index];
    segment.end = result.time[segment.end_index - 1];
    segment.series.resize(series_count);

    for (std::size_t j = 0; j < series_count; ++j) {
      SeriesSegment &entry = segment.series[j];
      entry.mean = mean_range(result.mean_series[j], segment.begin_index,
                              segment.end_index);

      std::vector<std::vector<double>> residuals(repeat_count,
                                                  std::vector<double>(point_count, 0.0));
      double pooled_square = 0.0;
      std::size_t pooled_degrees = 0;
      if (repeat_count >= 2) {
        for (std::size_t r = 0; r < repeat_count; ++r) {
          for (std::size_t i = segment.begin_index; i < segment.end_index; ++i) {
            residuals[r][i] = sets[r].series[j][i] - result.mean_series[j][i];
            pooled_square += residuals[r][i] * residuals[r][i];
          }
        }
        pooled_degrees = (repeat_count - 1) * segment.point_count;
        entry.noise_sigma = std::sqrt(pooled_square /
                                      std::max<std::size_t>(1, pooled_degrees));
      } else {
        const auto fit = linear_fit(result.time, result.mean_series[j],
                                    segment.begin_index, segment.end_index);
        entry.noise_sigma = residual_sigma(result.time, result.mean_series[j],
                                           segment.begin_index, segment.end_index,
                                           fit.first, fit.second);
        for (std::size_t i = segment.begin_index; i < segment.end_index; ++i) {
          residuals[0][i] = result.mean_series[j][i] - fit.first - fit.second * result.time[i];
        }
        pooled_degrees = segment.point_count > 2 ? segment.point_count - 2 : 1;
      }
      entry.rho1 = lag_one_correlation(residuals, segment.begin_index,
                                       segment.end_index);
      const double neff = effective_count(static_cast<double>(segment.point_count),
                                          entry.rho1);
      const double repeat_gain = std::sqrt(static_cast<double>(repeat_count));
      entry.mean_se = entry.noise_sigma / std::sqrt(neff) / repeat_gain;

      const int block_count = std::max(2, std::min(options.subblocks,
                                                   static_cast<int>(segment.point_count / 2)));
      std::vector<double> block_means;
      std::vector<double> block_ses;
      for (int block = 0; block < block_count; ++block) {
        const std::size_t begin = segment.begin_index +
            segment.point_count * static_cast<std::size_t>(block) / block_count;
        const std::size_t end = segment.begin_index +
            segment.point_count * static_cast<std::size_t>(block + 1) / block_count;
        block_means.push_back(mean_range(result.mean_series[j], begin, end));
        const double block_neff = effective_count(static_cast<double>(end - begin),
                                                  entry.rho1);
        block_ses.push_back(entry.noise_sigma / std::sqrt(block_neff) / repeat_gain);
      }
      const auto min_it = std::min_element(block_means.begin(), block_means.end());
      const auto max_it = std::max_element(block_means.begin(), block_means.end());
      const std::size_t min_index = static_cast<std::size_t>(min_it - block_means.begin());
      const std::size_t max_index = static_cast<std::size_t>(max_it - block_means.begin());
      const double amplitude = *max_it - *min_it;
      const double amplitude_se = std::sqrt(block_ses[min_index] * block_ses[min_index] +
                                            block_ses[max_index] * block_ses[max_index]);

      std::vector<double> repeat_slopes;
      for (std::size_t r = 0; r < repeat_count; ++r) {
        repeat_slopes.push_back(linear_fit(result.time, sets[r].series[j],
                                           segment.begin_index,
                                           segment.end_index).second);
      }
      const double slope = std::accumulate(repeat_slopes.begin(), repeat_slopes.end(), 0.0) /
                           static_cast<double>(repeat_count);
      double slope_se = 0.0;
      if (repeat_count >= 2) {
        double square = 0.0;
        for (double value : repeat_slopes) square += (value - slope) * (value - slope);
        slope_se = std::sqrt(square / static_cast<double>(repeat_count - 1)) /
                   repeat_gain;
      } else {
        const double mean_t = mean_range(result.time, segment.begin_index,
                                         segment.end_index);
        double sxx = 0.0;
        for (std::size_t i = segment.begin_index; i < segment.end_index; ++i) {
          const double centered = result.time[i] - mean_t;
          sxx += centered * centered;
        }
        slope_se = entry.noise_sigma / std::sqrt(std::max(kTiny, sxx));
        slope_se *= std::sqrt((1.0 + entry.rho1) /
                              std::max(kTiny, 1.0 - entry.rho1));
      }
      const double duration = std::max(kTiny, segment.end - segment.start);
      const int df = repeat_count >= 2 ? static_cast<int>(repeat_count - 1)
                                       : static_cast<int>(pooled_degrees);
      const double inflation = small_sample_inflation(df);
      entry.amplitude = make_metric(amplitude, amplitude_se, result.tolerances[j],
                                    inflation);
      entry.drift = make_metric(std::fabs(slope) * duration, slope_se * duration,
                                result.tolerances[j], inflation);
    }
  }

  for (std::size_t k = 0; k < result.segments.size(); ++k) {
    SegmentResult &segment = result.segments[k];
    segment.score = 100.0;
    for (std::size_t j = 0; j < series_count; ++j) {
      SeriesSegment &entry = segment.series[j];
      if (k + 1 < result.segments.size()) {
        const SeriesSegment &next = result.segments[k + 1].series[j];
        const double estimate = std::fabs(next.mean - entry.mean);
        const double standard_error = std::sqrt(entry.mean_se * entry.mean_se +
                                                next.mean_se * next.mean_se);
        const int df = repeat_count >= 2 ? static_cast<int>(repeat_count - 1)
                                         : static_cast<int>(segment.point_count - 2);
        entry.jump = make_metric(estimate, standard_error, result.tolerances[j],
                                 small_sample_inflation(df));
      } else {
        entry.jump = make_metric(0.0, kTiny, result.tolerances[j], 1.0);
      }
      entry.score = 100.0 * std::min(entry.amplitude.probability,
                          std::min(entry.drift.probability, entry.jump.probability));
      if (entry.score < segment.score) {
        segment.score = entry.score;
        segment.worst_series = static_cast<int>(j);
        const double p = std::min(entry.amplitude.probability,
                         std::min(entry.drift.probability, entry.jump.probability));
        if (p == entry.amplitude.probability) segment.limiting_metric = "amplitude";
        else if (p == entry.drift.probability) segment.limiting_metric = "drift";
        else segment.limiting_metric = "jump";
      }
    }
    segment.classification = classify(segment.score);
  }

  double running_tail = 100.0;
  for (std::size_t offset = 0; offset < result.segments.size(); ++offset) {
    const std::size_t k = result.segments.size() - 1 - offset;
    running_tail = std::min(running_tail, result.segments[k].score);
    result.segments[k].tail_score = running_tail;
  }

  const std::size_t eligible_count = result.segments.size() >=
      static_cast<std::size_t>(options.min_tail_segments)
      ? result.segments.size() - static_cast<std::size_t>(options.min_tail_segments) + 1
      : 0;
  result.best_eligible_tail_score = 0.0;
  for (std::size_t k = 0; k < eligible_count; ++k) {
    result.best_eligible_tail_score = std::max(result.best_eligible_tail_score,
                                               result.segments[k].tail_score);
    if (!std::isfinite(result.possible_start) && result.segments[k].tail_score >= 50.0) {
      result.possible_start = result.segments[k].start;
    }
    if (!std::isfinite(result.confirmed_start) &&
        result.segments[k].tail_score >= options.pass_probability * 100.0) {
      result.confirmed_start = result.segments[k].start;
    }
  }

  const double total_duration = std::max(kTiny, result.time.back() - result.time.front());
  double integral = 0.0;
  for (const SegmentResult &segment : result.segments) {
    integral += segment.tail_score * std::max(kTiny, segment.end - segment.start);
  }
  result.total_index = std::max(0.0, std::min(100.0, integral / total_duration));

  if (std::isfinite(result.confirmed_start)) result.status = "converged";
  else if (result.best_eligible_tail_score <= 5.0) result.status = "not_converged";
  else result.status = "inconclusive";

  if (repeat_count == 1) {
    result.warning = "Only one run was supplied; local residuals estimate noise and may mix noise with physical curvature.";
  } else if (repeat_count < 3) {
    result.warning = "Only two independent runs were supplied; small-sample inflation is applied and borderline results should receive a third run.";
  }
  if (!duplicate_paths.empty()) {
    if (!result.warning.empty()) result.warning += " ";
    result.warning += "Discarded " + std::to_string(duplicate_paths.size()) +
                      " byte-identical repeat(s); identical outputs are not independent Monte Carlo evidence.";
  }
  if (result.segments.size() >= 2) {
    double first_noise = 0.0;
    double last_noise = 0.0;
    for (const SeriesSegment &entry : result.segments.front().series) {
      first_noise += entry.noise_sigma * entry.noise_sigma;
    }
    for (const SeriesSegment &entry : result.segments.back().series) {
      last_noise += entry.noise_sigma * entry.noise_sigma;
    }
    first_noise = std::sqrt(first_noise / series_count);
    last_noise = std::sqrt(last_noise / series_count);
    if (first_noise > kTiny) result.noise_growth = last_noise / first_noise;
  }
  return result;
}

std::map<std::string, std::string> read_metadata(
    const std::vector<std::string> &files) {
  std::map<std::string, std::string> parameters;
  for (const std::string &file : files) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot open metadata file: " + file);
    std::string token;
    while (input >> token) {
      if (!token.empty() && token[0] == '#') token.erase(token.begin());
      const std::size_t equals = token.find('=');
      if (equals == std::string::npos || equals == 0 || equals + 1 >= token.size()) continue;
      std::string key = token.substr(0, equals);
      std::string value = token.substr(equals + 1);
      while (!value.empty() && (value.back() == ',' || value.back() == ';')) value.pop_back();
      parameters[key] = value;
    }
  }
  return parameters;
}

std::vector<ManifestCase> read_manifest(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open manifest: " + path);
  std::string line;
  if (!std::getline(input, line)) throw std::runtime_error("manifest is empty: " + path);
  const std::vector<std::string> header = parse_csv_row(line);
  int case_column = -1;
  int data_column = -1;
  std::vector<std::pair<int, std::string>> parameter_columns;
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == "case_id") case_column = static_cast<int>(i);
    else if (header[i] == "data_file") data_column = static_cast<int>(i);
    else if (header[i].find("param.") == 0) {
      parameter_columns.push_back({static_cast<int>(i), header[i].substr(6)});
    } else if (header[i].find("param:") == 0) {
      parameter_columns.push_back({static_cast<int>(i), header[i].substr(6)});
    }
  }
  if (case_column < 0 || data_column < 0) {
    throw std::runtime_error("manifest requires case_id and data_file columns");
  }
  std::map<std::string, ManifestCase> cases;
  const std::string base = parent_path(path);
  std::size_t line_number = 1;
  while (std::getline(input, line)) {
    ++line_number;
    if (trim(line).empty()) continue;
    const std::vector<std::string> cells = parse_csv_row(line);
    if (cells.size() != header.size()) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) +
                               ": CSV column count mismatch");
    }
    const std::string id = cells[case_column];
    if (id.empty() || cells[data_column].empty()) {
      throw std::runtime_error(path + ":" + std::to_string(line_number) +
                               ": case_id and data_file cannot be empty");
    }
    ManifestCase &entry = cases[id];
    entry.id = id;
    entry.files.push_back(join_path(base, cells[data_column]));
    for (const auto &column : parameter_columns) {
      const std::string &value = cells[column.first];
      auto found = entry.parameters.find(column.second);
      if (found != entry.parameters.end() && found->second != value) {
        throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                 ": parameter changes within case " + id);
      }
      entry.parameters[column.second] = value;
    }
  }
  std::vector<ManifestCase> result;
  for (auto &entry : cases) result.push_back(std::move(entry.second));
  if (result.empty()) throw std::runtime_error("manifest contains no cases");
  return result;
}

std::vector<double> solve_linear(std::vector<std::vector<double>> matrix,
                                 std::vector<double> rhs) {
  const std::size_t n = rhs.size();
  for (std::size_t pivot = 0; pivot < n; ++pivot) {
    std::size_t best = pivot;
    for (std::size_t row = pivot + 1; row < n; ++row) {
      if (std::fabs(matrix[row][pivot]) > std::fabs(matrix[best][pivot])) best = row;
    }
    if (std::fabs(matrix[best][pivot]) < 1.0e-12) {
      throw std::runtime_error("parameter tendency matrix is singular");
    }
    std::swap(matrix[pivot], matrix[best]);
    std::swap(rhs[pivot], rhs[best]);
    const double divisor = matrix[pivot][pivot];
    for (std::size_t col = pivot; col < n; ++col) matrix[pivot][col] /= divisor;
    rhs[pivot] /= divisor;
    for (std::size_t row = 0; row < n; ++row) {
      if (row == pivot) continue;
      const double factor = matrix[row][pivot];
      for (std::size_t col = pivot; col < n; ++col) {
        matrix[row][col] -= factor * matrix[pivot][col];
      }
      rhs[row] -= factor * rhs[pivot];
    }
  }
  return rhs;
}

TendencyModel fit_tendency(const std::vector<AnalysisResult> &results,
                           const Options &options) {
  TendencyModel model;
  if (results.size() < 3) {
    model.note = "At least three parameter cases are required for a tendency model.";
    return model;
  }
  std::set<std::string> common;
  for (const auto &item : results.front().parameters) {
    try {
      parse_double(item.second, item.first);
      common.insert(item.first);
    } catch (...) {
    }
  }
  for (std::size_t i = 1; i < results.size(); ++i) {
    for (auto it = common.begin(); it != common.end();) {
      auto value = results[i].parameters.find(*it);
      bool numeric = value != results[i].parameters.end();
      if (numeric) {
        try { parse_double(value->second, *it); } catch (...) { numeric = false; }
      }
      if (!numeric) it = common.erase(it);
      else ++it;
    }
  }
  if (common.empty()) {
    model.note = "No numeric parameter is shared by every case.";
    return model;
  }
  const std::vector<std::string> names(common.begin(), common.end());
  const std::size_t m = results.size();
  const std::size_t p = names.size();
  std::vector<std::vector<double>> x(m, std::vector<double>(p, 0.0));
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < p; ++j) {
      x[i][j] = parse_double(results[i].parameters.at(names[j]), names[j]);
    }
  }
  model.rows.resize(p);
  for (std::size_t j = 0; j < p; ++j) {
    double mean = 0.0;
    for (std::size_t i = 0; i < m; ++i) mean += x[i][j];
    mean /= static_cast<double>(m);
    double square = 0.0;
    for (std::size_t i = 0; i < m; ++i) square += (x[i][j] - mean) * (x[i][j] - mean);
    const double scale = std::sqrt(square / std::max<std::size_t>(1, m - 1));
    model.rows[j].parameter = names[j];
    model.rows[j].mean = mean;
    model.rows[j].scale = scale;
    for (std::size_t i = 0; i < m; ++i) {
      x[i][j] = scale > kTiny ? (x[i][j] - mean) / scale : 0.0;
    }
  }
  const std::size_t dimension = p + 1;
  std::vector<std::vector<double>> normal(dimension,
                                          std::vector<double>(dimension, 0.0));
  std::vector<double> rhs(dimension, 0.0);
  for (std::size_t i = 0; i < m; ++i) {
    std::vector<double> row(dimension, 1.0);
    for (std::size_t j = 0; j < p; ++j) row[j + 1] = x[i][j];
    for (std::size_t a = 0; a < dimension; ++a) {
      rhs[a] += row[a] * results[i].total_index;
      for (std::size_t b = 0; b < dimension; ++b) normal[a][b] += row[a] * row[b];
    }
  }
  for (std::size_t j = 1; j < dimension; ++j) normal[j][j] += options.ridge_penalty;
  std::vector<double> beta;
  try {
    beta = solve_linear(normal, rhs);
  } catch (const std::exception &error) {
    model.note = error.what();
    return model;
  }
  model.intercept = beta[0];
  model.coefficients.assign(beta.begin() + 1, beta.end());
  double residual_square = 0.0;
  for (std::size_t i = 0; i < m; ++i) {
    double prediction = beta[0];
    for (std::size_t j = 0; j < p; ++j) prediction += beta[j + 1] * x[i][j];
    const double residual = results[i].total_index - prediction;
    residual_square += residual * residual;
  }
  model.rmse = std::sqrt(residual_square / std::max<std::size_t>(1, m - dimension));
  for (std::size_t j = 0; j < p; ++j) {
    model.rows[j].coefficient = beta[j + 1];
    if (beta[j + 1] > 1.0) model.rows[j].direction = "higher tends to improve convergence";
    else if (beta[j + 1] < -1.0) model.rows[j].direction = "higher tends to reduce convergence";
    else model.rows[j].direction = "weak or unresolved tendency";
  }
  model.fitted = true;
  model.reliable = m >= std::max<std::size_t>(8, 2 * p + 4);
  model.note = model.reliable
      ? "Exploratory ridge tendency model; validate predictions with new runs."
      : "Too few parameter cases for a reliable nonlinear model; coefficients are screening hints only.";
  if (!options.predict_parameters.empty()) {
    bool complete = true;
    double prediction = beta[0];
    for (std::size_t j = 0; j < p; ++j) {
      auto found = options.predict_parameters.find(names[j]);
      if (found == options.predict_parameters.end()) {
        complete = false;
        break;
      }
      const double z = model.rows[j].scale > kTiny
                           ? (found->second - model.rows[j].mean) / model.rows[j].scale
                           : 0.0;
      prediction += beta[j + 1] * z;
    }
    if (complete) {
      model.predicted_index = std::max(0.0, std::min(100.0, prediction));
      model.predicted_low = std::max(0.0, model.predicted_index - 1.645 * model.rmse);
      model.predicted_high = std::min(100.0, model.predicted_index + 1.645 * model.rmse);
    } else {
      model.note += " Prediction omitted because --predict does not contain every fitted parameter.";
    }
  }
  return model;
}

std::string format_number(double value, int precision = 6) {
  if (!std::isfinite(value)) return "NA";
  std::ostringstream out;
  out << std::setprecision(precision) << value;
  return out.str();
}

std::string status_chinese(const std::string &status) {
  if (status == "converged") return "已确认收敛";
  if (status == "not_converged") return "未收敛";
  return "证据不足";
}

std::string score_color(double score) {
  if (score >= 95.0) return "#16803a";
  if (score >= 80.0) return "#4f9d3a";
  if (score >= 50.0) return "#d79518";
  return "#c43b3b";
}

std::vector<std::size_t> plot_indices(std::size_t count, int maximum) {
  std::vector<std::size_t> indices;
  const std::size_t kept = std::min<std::size_t>(count, std::max(2, maximum));
  indices.reserve(kept);
  for (std::size_t i = 0; i < kept; ++i) {
    indices.push_back(i * (count - 1) / (kept - 1));
  }
  return indices;
}

std::string svg_signal(const AnalysisResult &result, int max_points) {
  const int width = 920, height = 360;
  const int left = 70, right = 20, top = 24, bottom = 52;
  double ymin = std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  for (const auto &series : result.mean_series) {
    const double baseline = series.front();
    for (double value : series) {
      ymin = std::min(ymin, value - baseline);
      ymax = std::max(ymax, value - baseline);
    }
  }
  if (!(ymax > ymin)) { ymin -= 1.0; ymax += 1.0; }
  const double margin = 0.08 * (ymax - ymin);
  ymin -= margin;
  ymax += margin;
  const double xmin = result.time.front(), xmax = result.time.back();
  const auto indices = plot_indices(result.time.size(), max_points);
  static const char *colors[] = {"#1565c0", "#d84315", "#2e7d32", "#6a1b9a",
      "#00838f", "#ad1457", "#5d4037", "#455a64", "#9e9d24", "#ef6c00"};
  auto sx = [&](double x) { return left + (x - xmin) / std::max(kTiny, xmax - xmin) *
                                           (width - left - right); };
  auto sy = [&](double y) { return top + (ymax - y) / std::max(kTiny, ymax - ymin) *
                                          (height - top - bottom); };
  std::ostringstream out;
  out << "<svg viewBox='0 0 " << width << " " << height
      << "' role='img' aria-label='mean trajectories'>";
  out << "<rect width='100%' height='100%' fill='white'/>";
  for (int grid = 0; grid <= 4; ++grid) {
    const double y = ymin + (ymax - ymin) * grid / 4.0;
    out << "<line x1='" << left << "' x2='" << width-right << "' y1='" << sy(y)
        << "' y2='" << sy(y) << "' stroke='#e6eaf0'/><text x='" << left-8
        << "' y='" << sy(y)+4 << "' text-anchor='end' font-size='11' fill='#556'>"
        << html_escape(format_number(y, 3)) << "</text>";
  }
  for (std::size_t j = 0; j < result.mean_series.size(); ++j) {
    out << "<polyline fill='none' stroke='" << colors[j % 10]
        << "' stroke-width='1.5' points='";
    const double baseline = result.mean_series[j].front();
    for (std::size_t i : indices) {
      out << sx(result.time[i]) << ',' << sy(result.mean_series[j][i] - baseline) << ' ';
    }
    out << "'/><text x='" << (left + 8 + (j % 5) * 155) << "' y='"
        << (height - 15 - (j / 5) * 16) << "' font-size='11' fill='"
        << colors[j % 10] << "'>" << html_escape(result.series_names[j]) << "</text>";
  }
  out << "<line x1='" << left << "' x2='" << width-right << "' y1='" << height-bottom
      << "' y2='" << height-bottom << "' stroke='#334'/><text x='" << width/2
      << "' y='" << height-2 << "' text-anchor='middle' font-size='12'>time</text>"
      << "<text transform='translate(15," << height/2
      << ") rotate(-90)' text-anchor='middle' font-size='12'>value - value(t0)</text></svg>";
  return out.str();
}

std::string svg_segment_scores(const AnalysisResult &result) {
  const int width = 920, height = 300;
  const int left = 60, right = 24, top = 25, bottom = 60;
  const double bar_width = static_cast<double>(width - left - right) /
                           std::max<std::size_t>(1, result.segments.size());
  std::ostringstream out;
  out << "<svg viewBox='0 0 " << width << ' ' << height
      << "' role='img' aria-label='segment convergence scores'><rect width='100%' height='100%' fill='white'/>";
  for (int score : {0, 50, 80, 95, 100}) {
    const double y = top + (100.0 - score) / 100.0 * (height - top - bottom);
    out << "<line x1='" << left << "' x2='" << width-right << "' y1='" << y
        << "' y2='" << y << "' stroke='" << (score == 95 ? "#1b7837" : "#e6eaf0")
        << "' stroke-dasharray='" << (score == 95 ? "5 4" : "0") << "'/>"
        << "<text x='" << left-8 << "' y='" << y+4
        << "' text-anchor='end' font-size='11'>" << score << "</text>";
  }
  for (std::size_t k = 0; k < result.segments.size(); ++k) {
    const SegmentResult &segment = result.segments[k];
    const double x = left + k * bar_width + 0.13 * bar_width;
    const double h = segment.score / 100.0 * (height - top - bottom);
    const double y = height - bottom - h;
    out << "<rect x='" << x << "' y='" << y << "' width='" << 0.32*bar_width
        << "' height='" << h << "' fill='" << score_color(segment.score)
        << "'><title>SCI " << format_number(segment.score, 4) << "</title></rect>";
    const double th = segment.tail_score / 100.0 * (height - top - bottom);
    out << "<rect x='" << x + 0.38*bar_width << "' y='" << height-bottom-th
        << "' width='" << 0.32*bar_width << "' height='" << th
        << "' fill='#315b8a'><title>tail score " << format_number(segment.tail_score, 4)
        << "</title></rect><text x='" << left + (k+0.5)*bar_width << "' y='"
        << height-bottom+18 << "' text-anchor='middle' font-size='10'>"
        << html_escape(format_number(segment.start, 3) + "-" + format_number(segment.end, 3))
        << "</text>";
  }
  out << "<text x='" << width/2 << "' y='" << height-8
      << "' text-anchor='middle' font-size='12'>green/orange/red: segment SCI; blue: sustained-tail score</text></svg>";
  return out.str();
}

void write_segments_csv(const std::vector<AnalysisResult> &results,
                        const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write: " + path);
  out << "case_id,segment,start,end,points,SCI,tail_score,class,worst_series,limiting_metric,"
         "amplitude,amplitude_se,drift,drift_se,jump,jump_se,noise_sigma,rho1,tolerance\n";
  for (const AnalysisResult &result : results) {
    for (std::size_t k = 0; k < result.segments.size(); ++k) {
      const SegmentResult &segment = result.segments[k];
      const int j = std::max(0, segment.worst_series);
      const SeriesSegment &entry = segment.series[static_cast<std::size_t>(j)];
      out << csv_escape(result.case_id) << ',' << k + 1 << ',' << std::setprecision(12)
          << segment.start << ',' << segment.end << ',' << segment.point_count << ','
          << segment.score << ',' << segment.tail_score << ',' << segment.classification << ','
          << csv_escape(result.series_names[static_cast<std::size_t>(j)]) << ','
          << segment.limiting_metric << ',' << entry.amplitude.estimate << ','
          << entry.amplitude.standard_error << ',' << entry.drift.estimate << ','
          << entry.drift.standard_error << ',' << entry.jump.estimate << ','
          << entry.jump.standard_error << ',' << entry.noise_sigma << ',' << entry.rho1
          << ',' << result.tolerances[static_cast<std::size_t>(j)] << '\n';
    }
  }
}

void write_summary_csv(const std::vector<AnalysisResult> &results,
                       const std::string &path) {
  std::set<std::string> parameter_names;
  for (const AnalysisResult &result : results) {
    for (const auto &item : result.parameters) parameter_names.insert(item.first);
  }
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write: " + path);
  out << "case_id,status,total_index,best_tail_score,confirmed_start,possible_start,"
         "converged_end,repeats,points,series,noise_growth";
  for (const std::string &name : parameter_names) out << ',' << csv_escape("param." + name);
  out << '\n';
  for (const AnalysisResult &result : results) {
    out << csv_escape(result.case_id) << ',' << result.status << ',' << std::setprecision(12)
        << result.total_index << ',' << result.best_eligible_tail_score << ',';
    if (std::isfinite(result.confirmed_start)) out << result.confirmed_start;
    out << ',';
    if (std::isfinite(result.possible_start)) out << result.possible_start;
    out << ',' << result.time.back() << ',' << result.data_files.size() << ','
        << result.time.size() << ',' << result.mean_series.size() << ',';
    if (std::isfinite(result.noise_growth)) out << result.noise_growth;
    for (const std::string &name : parameter_names) {
      auto found = result.parameters.find(name);
      out << ',' << (found == result.parameters.end() ? "" : csv_escape(found->second));
    }
    out << '\n';
  }
}

void write_json(const std::vector<AnalysisResult> &results,
                const TendencyModel &model, const Options &options,
                const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write: " + path);
  out << std::setprecision(12) << "{\n  \"method\": \"heteroscedastic segmented equivalence\",\n"
      << "  \"absolute_tolerance\": " << options.absolute_tolerance << ",\n"
      << "  \"relative_tolerance\": " << options.relative_tolerance << ",\n"
      << "  \"pass_probability\": " << options.pass_probability << ",\n"
      << "  \"cases\": [\n";
  for (std::size_t i = 0; i < results.size(); ++i) {
    const AnalysisResult &result = results[i];
    out << "    {\"case_id\": \"" << json_escape(result.case_id) << "\", "
        << "\"status\": \"" << result.status << "\", "
        << "\"total_index\": " << result.total_index << ", "
        << "\"best_tail_score\": " << result.best_eligible_tail_score << ", "
        << "\"confirmed_start\": ";
    if (std::isfinite(result.confirmed_start)) out << result.confirmed_start;
    else out << "null";
    out << ", \"possible_start\": ";
    if (std::isfinite(result.possible_start)) out << result.possible_start;
    else out << "null";
    out << ", \"end\": " << result.time.back() << ", \"noise_growth\": ";
    if (std::isfinite(result.noise_growth)) out << result.noise_growth;
    else out << "null";
    out << ", \"parameters\": {";
    std::size_t parameter_index = 0;
    for (const auto &item : result.parameters) {
      if (parameter_index++) out << ',';
      out << "\"" << json_escape(item.first) << "\": \""
          << json_escape(item.second) << "\"";
    }
    out << "}, \"segments\": [";
    for (std::size_t k = 0; k < result.segments.size(); ++k) {
      if (k) out << ',';
      const SegmentResult &segment = result.segments[k];
      out << "{\"start\":" << segment.start << ",\"end\":" << segment.end
          << ",\"SCI\":" << segment.score << ",\"tail_score\":"
          << segment.tail_score << ",\"class\":\"" << segment.classification
          << "\",\"worst_series\":\""
          << json_escape(result.series_names[static_cast<std::size_t>(std::max(0, segment.worst_series))])
          << "\",\"limiting_metric\":\"" << segment.limiting_metric << "\"}";
    }
    out << "]}" << (i + 1 == results.size() ? "\n" : ",\n");
  }
  out << "  ],\n  \"parameter_tendency\": {\"fitted\": "
      << (model.fitted ? "true" : "false") << ", \"reliable\": "
      << (model.reliable ? "true" : "false") << ", \"rmse\": ";
  if (model.fitted) out << model.rmse;
  else out << "null";
  out << ", \"note\": \"" << json_escape(model.note) << "\", \"coefficients\": [";
  for (std::size_t j = 0; j < model.rows.size(); ++j) {
    if (j) out << ',';
    out << "{\"parameter\":\"" << json_escape(model.rows[j].parameter)
        << "\",\"standardized_coefficient\":" << model.rows[j].coefficient
        << ",\"direction\":\"" << json_escape(model.rows[j].direction) << "\"}";
  }
  out << "]}}\n";
}

void write_html(const std::vector<AnalysisResult> &results,
                const TendencyModel &model, const Options &options,
                const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write: " + path);
  out << "<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'>"
         "<meta name='viewport' content='width=device-width,initial-scale=1'>"
         "<title>仿真数据收敛报告</title><style>"
         "body{margin:0;background:#f4f6f8;color:#18212b;font-family:system-ui,-apple-system,'Segoe UI','Microsoft YaHei',sans-serif}"
         ".wrap{max-width:1180px;margin:auto;padding:28px}.hero{background:linear-gradient(135deg,#102a43,#1f6f8b);color:white;padding:28px;border-radius:18px}"
         ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;margin:18px 0}.card{background:white;border-radius:14px;padding:18px;box-shadow:0 2px 12px #102a4312}"
         ".big{font-size:30px;font-weight:750}.muted{color:#607080}.pill{display:inline-block;padding:5px 10px;border-radius:999px;color:white;font-weight:700}"
         "table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:9px;border-bottom:1px solid #e4e8ed}th{background:#f7f9fb;position:sticky;top:0}"
         ".scroll{overflow:auto}.warn{border-left:5px solid #d79518;background:#fff8e8;padding:12px 15px;border-radius:8px}"
         "svg{width:100%;height:auto}.refs li{margin:7px 0}a{color:#1565c0}code{background:#eef2f6;padding:2px 5px;border-radius:4px}"
         "</style></head><body><main class='wrap'><section class='hero'><h1>仿真数据收敛报告</h1>"
      << "<p>异方差分段等价检验；绝对容差 <code>" << options.absolute_tolerance
      << "</code>，相对容差 <code>" << options.relative_tolerance
      << "</code>，确认阈值 <code>" << options.pass_probability * 100.0
      << "%</code>。</p></section>";

  for (const AnalysisResult &result : results) {
    const std::string color = result.status == "converged" ? "#16803a" :
                              result.status == "not_converged" ? "#c43b3b" : "#d79518";
    out << "<section class='card'><h2>" << html_escape(result.case_id) << "</h2>"
        << "<span class='pill' style='background:" << color << "'>"
        << status_chinese(result.status) << "</span><div class='grid'>"
        << "<div><div class='muted'>总收敛指数 TCI</div><div class='big'>"
        << format_number(result.total_index, 4) << "</div></div>"
        << "<div><div class='muted'>最佳持续尾段分数</div><div class='big'>"
        << format_number(result.best_eligible_tail_score, 4) << "</div></div>"
        << "<div><div class='muted'>确认收敛区间</div><div class='big'>";
    if (std::isfinite(result.confirmed_start)) {
      out << '[' << format_number(result.confirmed_start, 5) << ", "
          << format_number(result.time.back(), 5) << ']';
    } else {
      out << "未确认";
    }
    out << "</div></div><div><div class='muted'>后段/前段噪声</div><div class='big'>"
        << (std::isfinite(result.noise_growth) ? format_number(result.noise_growth, 4) + "×" : "NA")
        << "</div></div></div>";
    if (!result.warning.empty()) out << "<p class='warn'>" << html_escape(result.warning) << "</p>";
    out << "<h3>初始参数</h3><div class='scroll'><table><tbody>";
    for (const auto &item : result.parameters) {
      out << "<tr><th>" << html_escape(item.first) << "</th><td>"
          << html_escape(item.second) << "</td></tr>";
    }
    out << "</tbody></table></div><h3>相对初值的潜在均值轨迹</h3>"
        << svg_signal(result, options.max_plot_points)
        << "<h3>分段 SCI 与持续尾段分数</h3>" << svg_segment_scores(result)
        << "<div class='scroll'><table><thead><tr><th>区间</th><th>SCI</th><th>持续尾段</th>"
           "<th>判断</th><th>最差序列</th><th>限制项</th></tr></thead><tbody>";
    for (const SegmentResult &segment : result.segments) {
      out << "<tr><td>" << format_number(segment.start, 4) << "–"
          << format_number(segment.end, 4) << "</td><td style='color:"
          << score_color(segment.score) << ";font-weight:700'>"
          << format_number(segment.score, 4) << "</td><td>"
          << format_number(segment.tail_score, 4) << "</td><td>"
          << html_escape(segment.classification) << "</td><td>"
          << html_escape(result.series_names[static_cast<std::size_t>(std::max(0, segment.worst_series))])
          << "</td><td>" << html_escape(segment.limiting_metric) << "</td></tr>";
    }
    out << "</tbody></table></div></section>";
  }

  if (model.fitted) {
    out << "<section class='card'><h2>初始参数—收敛倾向</h2><p class='warn'>"
        << html_escape(model.note) << "</p><p>系数表示参数增加一个样本标准差时，TCI 预计变化多少分；"
           "这是快速筛选模型，不代替新的独立验证。</p><table><thead><tr><th>参数</th><th>标准化系数</th><th>倾向</th></tr></thead><tbody>";
    for (const TendencyRow &row : model.rows) {
      out << "<tr><td>" << html_escape(row.parameter) << "</td><td>"
          << format_number(row.coefficient, 5) << "</td><td>"
          << html_escape(row.direction) << "</td></tr>";
    }
    out << "</tbody></table><p>训练残差 RMSE：" << format_number(model.rmse, 5) << " TCI 分。</p>";
    if (std::isfinite(model.predicted_index)) {
      out << "<p><b>输入参数的预测：</b>TCI=" << format_number(model.predicted_index, 5)
          << "，探索性 90% 区间 [" << format_number(model.predicted_low, 5) << ", "
          << format_number(model.predicted_high, 5) << "]。</p>";
    }
    out << "</section>";
  }

  out << "<section class='card'><h2>方法说明与参考模型</h2>"
         "<p>程序检验的是潜在均值在给定容差内是否持续平坦，而不是数据是否接近某个未知真值。"
         "重复运行用于分离真实趋势和随时间变化的噪声；单次运行则用稳健局部残差估计噪声，结论更保守。"
         "分段 SCI 是振幅、趋势和相邻段跳变三个等价概率的最小值；TCI 是持续尾段分数沿时间的归一化积分。</p>"
         "<ol class='refs'>"
         "<li>Heidelberger & Welch (1983), <a href='https://doi.org/10.1287/opre.31.6.1109'>Simulation Run Length Control in the Presence of an Initial Transient</a>.</li>"
         "<li>Flegal & Jones (2010), <a href='https://doi.org/10.1214/09-AOS735'>Batch Means and Spectral Variance Estimators in MCMC</a>.</li>"
         "<li>Górecki, Horváth & Kokoszka (2018), <a href='https://doi.org/10.1016/j.ecosta.2017.07.005'>Change Point Detection in Heteroscedastic Time Series</a>.</li>"
         "<li>Kersting et al. (2007), <a href='https://doi.org/10.1145/1273496.1273546'>Most Likely Heteroscedastic Gaussian Process Regression</a>.</li>"
         "<li>Glynn & Whitt (1992), <a href='https://web.stanford.edu/~glynn/papers/1992/GW92a.html'>Sequential Stopping Rules for Stochastic Simulations</a>.</li>"
         "</ol></section></main></body></html>";
}

void print_result(const AnalysisResult &result) {
  std::cout << "case=" << result.case_id << "\n"
            << "status=" << result.status << " (" << status_chinese(result.status) << ")\n"
            << "TCI=" << format_number(result.total_index, 6) << "\n"
            << "best_tail_score=" << format_number(result.best_eligible_tail_score, 6) << "\n"
            << "convergence_interval=";
  if (std::isfinite(result.confirmed_start)) {
    std::cout << '[' << result.confirmed_start << ',' << result.time.back() << "]\n";
  } else {
    std::cout << "not_confirmed_by_t=" << result.time.back() << "\n";
  }
  if (std::isfinite(result.possible_start)) {
    std::cout << "possible_onset=" << result.possible_start << "\n";
  }
  if (std::isfinite(result.noise_growth)) {
    std::cout << "noise_growth_last_vs_first=" << result.noise_growth << "x\n";
  }
  if (!result.warning.empty()) std::cout << "warning=" << result.warning << "\n";
}

void usage(std::ostream &out) {
  out << "Fast convergence analyzer for heteroscedastic simulation data\n\n"
      << "Single case:\n"
      << "  convergence_analyzer --data run1.dat --data run2.dat [options]\n\n"
      << "Batch/parameter tendency:\n"
      << "  convergence_analyzer --manifest cases.csv [--predict name=value ...] [options]\n\n"
      << "Core options:\n"
      << "  --param NAME=VALUE       Initial parameter; repeat for n dimensions\n"
      << "  --metadata FILE          Read key=value tokens from program.out/run-info.txt\n"
      << "  --time-col N             Zero-based time column (default 0)\n"
      << "  --value-cols SPEC        Zero-based columns, e.g. 1 or 4:13 or 4: (default 4:)\n"
      << "  --series-names A,B       Labels for selected columns\n"
      << "  --segments N             Equal-duration segments (default 5)\n"
      << "  --segment-width X        Physical segment width; overrides --segments\n"
      << "  --atol X                 Absolute practical tolerance (default 1e-4)\n"
      << "  --rtol X                 Fraction of robust full-curve range (default 0.01)\n"
      << "  --pass-prob X            Confirmation probability (default 0.95)\n"
      << "  --min-tail-segments N    Required consecutive tail segments (default 2)\n"
      << "  --out-prefix PATH        Output prefix (default convergence)\n"
      << "  --self-test              Run built-in synthetic tests\n\n"
      << "Manifest columns: case_id,data_file,param.wc,param.eta,...\n";
}

Options parse_options(int argc, char **argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    auto value = [&](const std::string &name) -> std::string {
      if (i + 1 >= argc) throw std::runtime_error(name + " requires a value");
      return argv[++i];
    };
    if (argument == "--data") options.data_files.push_back(value(argument));
    else if (argument == "--manifest") options.manifest = value(argument);
    else if (argument == "--metadata") options.metadata_files.push_back(value(argument));
    else if (argument == "--param") {
      const auto item = parse_key_value(value(argument), argument);
      options.parameters[item.first] = item.second;
    } else if (argument == "--predict") {
      const auto item = parse_key_value(value(argument), argument);
      options.predict_parameters[item.first] = parse_double(item.second, item.first);
    } else if (argument == "--time-col") options.time_column = parse_int(value(argument), argument);
    else if (argument == "--value-cols") options.value_columns = value(argument);
    else if (argument == "--series-names") options.series_names = split(value(argument), ',');
    else if (argument == "--segments") options.segment_count = parse_int(value(argument), argument);
    else if (argument == "--segment-width") options.segment_width = parse_double(value(argument), argument);
    else if (argument == "--subblocks") options.subblocks = parse_int(value(argument), argument);
    else if (argument == "--min-tail-segments") options.min_tail_segments = parse_int(value(argument), argument);
    else if (argument == "--max-plot-points") options.max_plot_points = parse_int(value(argument), argument);
    else if (argument == "--atol") options.absolute_tolerance = parse_double(value(argument), argument);
    else if (argument == "--rtol") options.relative_tolerance = parse_double(value(argument), argument);
    else if (argument == "--pass-prob") options.pass_probability = parse_double(value(argument), argument);
    else if (argument == "--ridge") options.ridge_penalty = parse_double(value(argument), argument);
    else if (argument == "--out-prefix") options.out_prefix = value(argument);
    else if (argument == "--self-test") options.self_test = true;
    else if (argument == "--help" || argument == "-h") {
      usage(std::cout);
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (options.segment_count < 1 || options.subblocks < 2 ||
      options.min_tail_segments < 1 || options.max_plot_points < 2) {
    throw std::runtime_error("segment and plotting counts must be positive");
  }
  if (options.absolute_tolerance < 0.0 || options.relative_tolerance < 0.0 ||
      options.absolute_tolerance + options.relative_tolerance <= 0.0) {
    throw std::runtime_error("at least one positive tolerance is required");
  }
  if (!(options.pass_probability > 0.5 && options.pass_probability < 1.0)) {
    throw std::runtime_error("--pass-prob must be between 0.5 and 1");
  }
  return options;
}

void run_self_test() {
  const std::size_t n = 501;
  DataSet converged, drifting;
  converged.path = "synthetic_converged";
  drifting.path = "synthetic_drifting";
  converged.series.assign(1, std::vector<double>(n));
  drifting.series.assign(1, std::vector<double>(n));
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / 5.0;
    converged.time.push_back(t);
    drifting.time.push_back(t);
    const double deterministic_noise = 0.0015 * std::sin(1.7 * t) +
                                       0.0008 * std::sin(3.1 * t + 0.4);
    converged.series[0][i] = 1.0 + std::exp(-t / 12.0) + deterministic_noise;
    drifting.series[0][i] = 1.0 + 0.002 * t + deterministic_noise;
  }
  const double tail_converged = std::fabs(converged.series[0].back() -
                                          converged.series[0][400]);
  const double tail_drifting = std::fabs(drifting.series[0].back() -
                                         drifting.series[0][400]);
  if (!(tail_converged < 0.01 && tail_drifting > 0.03)) {
    throw std::runtime_error("built-in synthetic test failed");
  }
  std::cout << "self-test passed: converged tail change=" << tail_converged
            << ", drifting tail change=" << tail_drifting << '\n';
}

}  // namespace

int main(int argc, char **argv) {
  try {
    Options options = parse_options(argc, argv);
    if (options.self_test) {
      run_self_test();
      return 0;
    }
    if (!options.manifest.empty() && !options.data_files.empty()) {
      throw std::runtime_error("use either --manifest or --data, not both");
    }
    if (options.manifest.empty() && options.data_files.empty()) {
      usage(std::cerr);
      throw std::runtime_error("no input data supplied");
    }

    std::vector<AnalysisResult> results;
    if (!options.manifest.empty()) {
      const std::vector<ManifestCase> cases = read_manifest(options.manifest);
      for (const ManifestCase &entry : cases) {
        results.push_back(analyze_case(entry.id, entry.files, entry.parameters, options));
      }
    } else {
      std::map<std::string, std::string> parameters = read_metadata(options.metadata_files);
      for (const auto &item : options.parameters) parameters[item.first] = item.second;
      results.push_back(analyze_case("single_case", options.data_files, parameters, options));
    }

    const TendencyModel tendency = fit_tendency(results, options);
    write_summary_csv(results, options.out_prefix + "_summary.csv");
    write_segments_csv(results, options.out_prefix + "_segments.csv");
    write_json(results, tendency, options, options.out_prefix + "_summary.json");
    write_html(results, tendency, options, options.out_prefix + "_report.html");
    for (const AnalysisResult &result : results) print_result(result);
    if (tendency.fitted) {
      std::cout << "parameter_model=" << tendency.note << "\n";
      if (std::isfinite(tendency.predicted_index)) {
        std::cout << "predicted_TCI=" << tendency.predicted_index << " ["
                  << tendency.predicted_low << ',' << tendency.predicted_high << "]\n";
      }
    }
    std::cout << "outputs=" << options.out_prefix << "_report.html, "
              << options.out_prefix << "_summary.csv, " << options.out_prefix
              << "_segments.csv, " << options.out_prefix << "_summary.json\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 2;
  }
}

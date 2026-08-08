# 快速仿真数据收敛分析器

`convergence_analyzer.cpp` 是一个无第三方运行时依赖的 C++17 命令行工具，适合本项目的
SepMB 输出，也可用于普通的“一列时间 + 一列结果”数据。

它输出：

- 每个时间段的收敛指数 `SCI`；
- 从该时间段直到结束的持续尾段分数；
- 总收敛指数 `TCI`；
- 确认收敛区间与可能的收敛起点；
- 后段相对前段的噪声增长倍数；
- 多组初始参数的收敛倾向；
- CSV 和 JSON 数值摘要；仅在显式传入 `--html-report` 时生成 HTML。

## 编译

Linux/云服务器：

```bash
make -f Makefile.convergence
```

Windows MinGW：

```powershell
g++ -O3 -std=c++17 -DNDEBUG convergence_analyzer.cpp -o convergence_analyzer.exe
```

程序只使用 C++ 标准库，数据分析为线性扫描。默认流程不生成网页或图片。

也可以直接使用启动脚本；它们会在源码变化后自动重新编译：

```powershell
.\run_convergence.cmd --help
```

```bash
bash run_convergence.sh --help
```

## 直接分析本项目的一组参数

SepMB 文件中第 0 列是时间，第 4 列开始是各轨道占据数：

```powershell
.\convergence_analyzer.exe `
  --data scan_cloud_t50_20260808_grid_optimized\wc6.5_eta7.5e-5_rep1\ahm-sepmb-s10-n5-2000000.dat `
  --data scan_cloud_t50_20260808_grid_optimized\wc6.5_eta7.5e-5_rep2\ahm-sepmb-s10-n5-2000000.dat `
  --data scan_cloud_t50_20260808_grid_optimized\wc6.5_eta7.5e-5_rep3\ahm-sepmb-s10-n5-2000000.dat `
  --param wc_eV=6.5 --param eta=7.5e-5 --param ntraj=2000000 `
  --time-col 0 --value-cols 4:13 --segment-width 10 `
  --atol 1e-4 --rtol 0 --out-prefix tmp\wc65_convergence
```

若只分析一个一维结果列，例如普通文件中的第 1 列：

```powershell
.\convergence_analyzer.exe --data result.dat --time-col 0 --value-cols 1 `
  --param alpha=0.2 --param beta=4.0 --atol 0.001 --out-prefix result_check
```

`--metadata program.out` 或 `--metadata run-info.txt` 可以自动读取其中的 `key=value`
初始参数；命令行 `--param` 会覆盖同名值。

## 批量参数分析

清单格式：

```csv
case_id,data_file,param.wc_eV,param.eta,param.ntraj
case_a,case_a/run1.dat,6.5,7.5e-5,2000000
case_a,case_a/run2.dat,6.5,7.5e-5,2000000
case_b,case_b/run1.dat,4.6,1.3e-4,2000000
case_b,case_b/run2.dat,4.6,1.3e-4,2000000
```

相同 `case_id` 的行视为独立重复，所有 `param.` 开头的列都作为任意维初始参数：

```powershell
.\convergence_analyzer.exe --manifest convergence_manifest.project.csv `
  --time-col 0 --value-cols 4:13 --segment-width 10 `
  --atol 1e-4 --rtol 0 --out-prefix tmp\project_convergence
```

还可以输入一个新参数点，得到探索性的 TCI 预测区间：

```powershell
.\convergence_analyzer.exe --manifest convergence_manifest.project.csv `
  --predict wc_eV=5.0 --predict eta=1.0e-4 --predict ntraj=2000000 `
  --value-cols 4:13 --segment-width 10 --atol 1e-4
```

当参数组合数少于 `max(8, 2n+4)` 时，报告会明确标记参数倾向模型不可靠，不会把少量
样本的回归结果冒充确定结论。

## 振荡信号的参数筛选

`SCI/TCI` 检验的是时间序列是否趋于平坦，不适合直接淘汰本项目中应当保留的量子振荡。
参数筛选使用独立重复之间的分段信噪比 `Q`：

```powershell
python rank_sci_q.py --manifest manifest.csv --out-prefix ranking `
  --segment-width 10 --required-until 60 --min-q 2 --min-repeats 3
```

该命令只输出 CSV/JSON。所有完整分段都满足 `Q >= 2`、独立重复数不少于 3 且电子数
守恒误差不超过 `1e-8` 时，所要求的时间区间才标记为有效。

绘图由独立模块按需完成：

```powershell
python plot_scan_results.py --manifest manifest.csv --ranking ranking.csv `
  --out-prefix result --top 5
```

## 指标含义

每段分别检验：

1. 子块均值的最大跨度；
2. 线性趋势在该段造成的总变化；
3. 与下一段均值的跳变。

三项都小于容差的保守等价概率取最小值，得到 `SCI`。程序从后向前取 SCI 的最小值，
得到持续尾段分数，防止中间偶然平静的区间被误认成收敛。`TCI` 是持续尾段分数沿时间
的归一化积分，范围为 0–100。

- `SCI >= 95`：强收敛；
- `80 <= SCI < 95`：可用但未达到默认确认线；
- `50 <= SCI < 80`：证据不足；
- `SCI < 50`：未收敛。

至少需要连续两个尾段超过确认阈值，才输出确认收敛区间。只有一次运行时，噪声由局部
稳健残差估计，会混入真实曲率；若程序输出“证据不足”，优先增加一次独立重复。

## 方法参考

- Heidelberger & Welch (1983), *Simulation Run Length Control in the Presence of an Initial Transient*, Operations Research. DOI: 10.1287/opre.31.6.1109.
- Flegal & Jones (2010), *Batch Means and Spectral Variance Estimators in Markov Chain Monte Carlo*, Annals of Statistics. DOI: 10.1214/09-AOS735.
- Górecki, Horváth & Kokoszka (2018), *Change Point Detection in Heteroscedastic Time Series*, Econometrics and Statistics. DOI: 10.1016/j.ecosta.2017.07.005.
- Kersting et al. (2007), *Most Likely Heteroscedastic Gaussian Process Regression*, ICML. DOI: 10.1145/1273496.1273546.
- Glynn & Whitt (1992), *The Asymptotic Validity of Sequential Stopping Rules for Stochastic Simulations*, Annals of Applied Probability.

本程序为了速度和无依赖部署，采用上述理论启发的分段等价近似，并没有声称实现完整的
异方差 GP、PELT 或贝叶斯生存模型。参数趋势使用带岭惩罚的标准化线性筛选模型；正式
外推必须用新的独立仿真验证。

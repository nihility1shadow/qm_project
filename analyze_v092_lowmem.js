#!/usr/bin/env node
"use strict";

const fs = require("fs");
const path = require("path");

function argument(name, fallback = undefined) {
  const index = process.argv.indexOf(`--${name}`);
  return index >= 0 ? process.argv[index + 1] : fallback;
}

function loadData(file) {
  return fs.readFileSync(file, "utf8")
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line && !line.startsWith("#"))
    .map((line) => line.split(/\s+/).map(Number));
}

function interpolate(rows, time, firstOccupationColumn, norb, cursor) {
  while (cursor.value + 1 < rows.length && rows[cursor.value + 1][0] < time) {
    cursor.value += 1;
  }
  const left = rows[cursor.value];
  const right = rows[Math.min(cursor.value + 1, rows.length - 1)];
  if (right[0] === left[0]) return left.slice(firstOccupationColumn, firstOccupationColumn + norb);
  const alpha = (time - left[0]) / (right[0] - left[0]);
  return Array.from({ length: norb }, (_, orbital) =>
    left[firstOccupationColumn + orbital] * (1 - alpha)
      + right[firstOccupationColumn + orbital] * alpha);
}

function rms(values) {
  return Math.sqrt(values.reduce((sum, value) => sum + value * value, 0) / values.length);
}

function main() {
  const qmFile = argument("qm");
  const simulationFile = argument("simulation");
  const outputDir = argument("out-dir");
  const tmax = Number(argument("tmax"));
  const windowWidth = Number(argument("window", "25"));
  const ntraj = Number(argument("ntraj"));
  const targetQ = Number(argument("target-q", "10"));
  const norb = Number(argument("norb", "10"));
  const active = argument("active", "0,5,6,7,8,9").split(",").map(Number);
  if (!qmFile || !simulationFile || !outputDir || !Number.isFinite(tmax)) {
    throw new Error("Need --qm, --simulation, --out-dir, and --tmax");
  }

  const qm = loadData(qmFile);
  const simulation = loadData(simulationFile).filter((row) => row[0] <= tmax + 1e-12);
  const initial = qm[0].slice(4, 4 + norb);
  const cursor = { value: 0 };
  const samples = simulation.map((row) => {
    const reference = interpolate(qm, row[0], 4, norb, cursor);
    const estimate = row.slice(4, 4 + norb);
    return {
      time: row[0],
      signal: reference.map((value, orbital) => value - initial[orbital]),
      error: estimate.map((value, orbital) => value - reference[orbital]),
    };
  });

  const windows = [];
  for (let start = 0; start < tmax - 1e-12; start += windowWidth) {
    const stop = Math.min(start + windowWidth, tmax);
    const selected = samples.filter((sample) =>
      sample.time >= start - 1e-12 && sample.time <= stop + 1e-12);
    const signal = [];
    const error = [];
    const orbitalQ = {};
    for (const orbital of active) {
      const orbitalSignal = selected.map((sample) => sample.signal[orbital]);
      const orbitalError = selected.map((sample) => sample.error[orbital]);
      orbitalQ[orbital] = rms(orbitalSignal) / Math.max(rms(orbitalError), 1e-300);
      signal.push(...orbitalSignal);
      error.push(...orbitalError);
    }
    const aggregateQ = rms(signal) / Math.max(rms(error), 1e-300);
    const minOrbitalQ = Math.min(...Object.values(orbitalQ));
    windows.push({
      start,
      stop,
      points: selected.length,
      aggregate_Q: aggregateQ,
      min_orbital_Q: minOrbitalQ,
      orbital_Q: orbitalQ,
      projected_ntraj_aggregate_Q10: Number.isFinite(ntraj)
        ? ntraj * (targetQ / aggregateQ) ** 2 : null,
      projected_ntraj_strict_Q10: Number.isFinite(ntraj)
        ? ntraj * (targetQ / minOrbitalQ) ** 2 : null,
    });
  }

  const summary = {
    qm: path.resolve(qmFile),
    simulation: path.resolve(simulationFile),
    ntraj,
    tmax: samples.at(-1).time,
    active_orbitals: active,
    minimum_aggregate_Q: Math.min(...windows.map((row) => row.aggregate_Q)),
    minimum_strict_Q: Math.min(...windows.map((row) => row.min_orbital_Q)),
    windows,
  };
  fs.mkdirSync(outputDir, { recursive: true });
  fs.writeFileSync(path.join(outputDir, "q_metrics.json"), `${JSON.stringify(summary, null, 2)}\n`);
  const csv = [
    "start,stop,points,aggregate_Q,min_orbital_Q,projected_ntraj_aggregate_Q10,projected_ntraj_strict_Q10",
    ...windows.map((row) => [row.start, row.stop, row.points, row.aggregate_Q,
      row.min_orbital_Q, row.projected_ntraj_aggregate_Q10,
      row.projected_ntraj_strict_Q10].join(",")),
  ].join("\n");
  fs.writeFileSync(path.join(outputDir, "q_windows.csv"), `${csv}\n`);
  for (const row of windows) {
    console.log(`${row.start}-${row.stop}: aggregate_Q=${row.aggregate_Q.toPrecision(6)} `
      + `strict_Q=${row.min_orbital_Q.toPrecision(6)} `
      + `strict_N10=${Math.ceil(row.projected_ntraj_strict_Q10)}`);
  }
}

main();

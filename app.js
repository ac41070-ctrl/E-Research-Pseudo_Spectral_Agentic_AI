const N = 1024;
const L = 100;
const Y_MIN = -0.4;
const Y_MAX = 1.2;
const dx = (2 * L) / N;
const x = Array.from({ length: N }, (_, i) => -L + i * dx);
const k = Array.from({ length: N }, (_, i) => {
  const n = i < N / 2 ? i : i - N;
  return (Math.PI * n) / L;
});

const state = {
  mode: "harmonic",
  running: true,
  t: 0,
  psiRe: new Float64Array(N),
  psiIm: new Float64Array(N),
  potential: new Float64Array(N),
};

const els = {
  waveCanvas: document.getElementById("waveCanvas"),
  probabilityCanvas: document.getElementById("probabilityCanvas"),
  energyPositionCanvas: document.getElementById("energyPositionCanvas"),
  playBtn: document.getElementById("playBtn"),
  resetBtn: document.getElementById("resetBtn"),
  scenarioText: document.getElementById("scenarioText"),
  timeReadout: document.getElementById("timeReadout"),
  normReadout: document.getElementById("normReadout"),
  energyReadout: document.getElementById("energyReadout"),
  meanXReadout: document.getElementById("meanXReadout"),
  meanPReadout: document.getElementById("meanPReadout"),
  widthReadout: document.getElementById("widthReadout"),
  alertReadout: document.getElementById("alertReadout"),
  stripTime: document.getElementById("stripTime"),
  stripNorm: document.getElementById("stripNorm"),
  stripEnergy: document.getElementById("stripEnergy"),
  stripMeanX: document.getElementById("stripMeanX"),
  stripMeanP: document.getElementById("stripMeanP"),
  stripWidth: document.getElementById("stripWidth"),
  omegaRow: document.getElementById("omegaRow"),
  barrierHeightRow: document.getElementById("barrierHeightRow"),
  barrierWidthRow: document.getElementById("barrierWidthRow"),
};

const controls = Object.fromEntries(
  ["x0", "k0", "sigma", "omega", "barrierHeight", "barrierWidth", "dt", "speed"].map((id) => [
    id,
    document.getElementById(id),
  ]),
);

function fft(re, im, inverse = false) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i += 1) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
  }

  for (let len = 2; len <= n; len <<= 1) {
    const angle = ((inverse ? 2 : -2) * Math.PI) / len;
    const wLenRe = Math.cos(angle);
    const wLenIm = Math.sin(angle);
    for (let i = 0; i < n; i += len) {
      let wRe = 1;
      let wIm = 0;
      for (let j = 0; j < len / 2; j += 1) {
        const uRe = re[i + j];
        const uIm = im[i + j];
        const vRe = re[i + j + len / 2] * wRe - im[i + j + len / 2] * wIm;
        const vIm = re[i + j + len / 2] * wIm + im[i + j + len / 2] * wRe;
        re[i + j] = uRe + vRe;
        im[i + j] = uIm + vIm;
        re[i + j + len / 2] = uRe - vRe;
        im[i + j + len / 2] = uIm - vIm;
        const nextRe = wRe * wLenRe - wIm * wLenIm;
        wIm = wRe * wLenIm + wIm * wLenRe;
        wRe = nextRe;
      }
    }
  }

  if (inverse) {
    for (let i = 0; i < n; i += 1) {
      re[i] /= n;
      im[i] /= n;
    }
  }
}

function value(id) {
  return Number(controls[id].value);
}

function setDefaultsForMode(mode) {
  if (mode === "harmonic") {
    controls.x0.value = -18;
    controls.k0.value = 0;
    controls.sigma.value = 2;
    controls.omega.value = 0.16;
  } else if (mode === "tunnel") {
    controls.x0.value = -35;
    controls.k0.value = 2.0;
    controls.sigma.value = 4;
    controls.barrierHeight.value = 2.3;
    controls.barrierWidth.value = 3;
  } else {
    controls.x0.value = -35;
    controls.k0.value = 1.2;
    controls.sigma.value = 3;
  }
  updateOutputs();
}

function updateOutputs() {
  Object.entries(controls).forEach(([id, el]) => {
    const out = document.getElementById(`${id}Out`);
    if (!out) return;
    out.textContent = id === "speed" ? `${el.value}x` : Number(el.value).toFixed(id === "dt" ? 3 : 2);
  });
}

function updatePotential() {
  const omega = value("omega");
  const height = value("barrierHeight");
  const width = value("barrierWidth");
  for (let i = 0; i < N; i += 1) {
    const xi = x[i];
    if (state.mode === "harmonic") state.potential[i] = 0.5 * omega * omega * xi * xi;
    else if (state.mode === "tunnel") state.potential[i] = Math.abs(xi) < width / 2 ? height : 0;
    else state.potential[i] = 0;
  }

  els.omegaRow.hidden = state.mode !== "harmonic";
  els.barrierHeightRow.hidden = state.mode !== "tunnel";
  els.barrierWidthRow.hidden = state.mode !== "tunnel";
  els.scenarioText.textContent = {
    harmonic: "A Gaussian packet oscillates in a quadratic trap. The center should move back and forth while norm and energy remain nearly constant.",
    tunnel: "A moving packet strikes a finite barrier. Reflected and transmitted components reveal tunneling and interference.",
    free: "A free packet spreads because different momentum components carry different phases.",
  }[state.mode];
}

function initializeWavepacket() {
  const x0 = value("x0");
  const k0 = value("k0");
  const sigma = value("sigma");
  const amp = (1 / (2 * Math.PI * sigma * sigma)) ** 0.25;
  for (let i = 0; i < N; i += 1) {
    const envelope = amp * Math.exp(-((x[i] - x0) ** 2) / (4 * sigma * sigma));
    const phase = k0 * (x[i] - x0);
    state.psiRe[i] = envelope * Math.cos(phase);
    state.psiIm[i] = envelope * Math.sin(phase);
  }
  normalize();
  state.t = 0;
}

function normalize() {
  const norm = Math.sqrt(measure().norm);
  if (!Number.isFinite(norm) || norm === 0) return;
  for (let i = 0; i < N; i += 1) {
    state.psiRe[i] /= norm;
    state.psiIm[i] /= norm;
  }
}

function complexPhaseMultiply(re, im, angle) {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [re * c - im * s, re * s + im * c];
}

function step() {
  const dt = value("dt");

  for (let i = 0; i < N; i += 1) {
    [state.psiRe[i], state.psiIm[i]] = complexPhaseMultiply(
      state.psiRe[i],
      state.psiIm[i],
      -0.5 * state.potential[i] * dt,
    );
  }

  fft(state.psiRe, state.psiIm, false);
  for (let i = 0; i < N; i += 1) {
    [state.psiRe[i], state.psiIm[i]] = complexPhaseMultiply(
      state.psiRe[i],
      state.psiIm[i],
      -0.5 * k[i] * k[i] * dt,
    );
  }
  fft(state.psiRe, state.psiIm, true);

  for (let i = 0; i < N; i += 1) {
    [state.psiRe[i], state.psiIm[i]] = complexPhaseMultiply(
      state.psiRe[i],
      state.psiIm[i],
      -0.5 * state.potential[i] * dt,
    );
  }

  state.t += dt;
}

function measure() {
  let norm = 0;
  let meanX = 0;
  let meanX2 = 0;
  let meanP = 0;
  let potentialEnergy = 0;
  for (let i = 0; i < N; i += 1) {
    const p = state.psiRe[i] ** 2 + state.psiIm[i] ** 2;
    const next = (i + 1) % N;
    const prev = (i - 1 + N) % N;
    const dReDx = (state.psiRe[next] - state.psiRe[prev]) / (2 * dx);
    const dImDx = (state.psiIm[next] - state.psiIm[prev]) / (2 * dx);
    norm += p * dx;
    meanX += x[i] * p * dx;
    meanX2 += x[i] * x[i] * p * dx;
    meanP += (state.psiRe[i] * dImDx - state.psiIm[i] * dReDx) * dx;
    potentialEnergy += state.potential[i] * p * dx;
  }

  const reK = Float64Array.from(state.psiRe);
  const imK = Float64Array.from(state.psiIm);
  fft(reK, imK, false);
  let kineticEnergy = 0;
  for (let i = 0; i < N; i += 1) {
    const scale = dx / N;
    const pk = (reK[i] ** 2 + imK[i] ** 2) * scale;
    kineticEnergy += 0.5 * k[i] * k[i] * pk;
  }

  return {
    norm,
    meanX: meanX / Math.max(norm, 1e-12),
    meanP: meanP / Math.max(norm, 1e-12),
    width: Math.sqrt(Math.max(meanX2 / Math.max(norm, 1e-12) - (meanX / Math.max(norm, 1e-12)) ** 2, 0)),
    energy: kineticEnergy + potentialEnergy,
  };
}

function drawWave() {
  const canvas = els.waveCanvas;
  const ctx = canvas.getContext("2d");
  const w = canvas.width;
  const h = canvas.height;
  ctx.clearRect(0, 0, w, h);
  drawGrid(ctx, w, h);
  drawWaveAxes(ctx, w, h);

  const modulus = Array.from({ length: N }, (_, i) => Math.hypot(state.psiRe[i], state.psiIm[i]));
  const maxV = potentialDisplayScale();

  drawPotentialOverlay(ctx, w, h, maxV);
  drawWaveCurve(ctx, Array.from(state.psiRe), "#74d3ae", w, h);
  drawWaveCurve(ctx, Array.from(state.psiIm), "#ff8f70", w, h);
  drawWaveCurve(ctx, modulus, "#7fb5ff", w, h, true);
  const measurement = measure();
  drawWidthMarker(ctx, measurement, w, h);
  drawWaveReadout(ctx, measurement, w);
}

function plotBounds(w, h) {
  return {
    left: 54,
    right: w - 22,
    top: 24,
    bottom: h - 42,
  };
}

function xToPixel(xValue, w, h) {
  const bounds = plotBounds(w, h);
  return bounds.left + ((xValue + L) / (2 * L)) * (bounds.right - bounds.left);
}

function yToPixel(yValue, w, h) {
  const bounds = plotBounds(w, h);
  return bounds.bottom - ((yValue - Y_MIN) / (Y_MAX - Y_MIN)) * (bounds.bottom - bounds.top);
}

function drawWaveAxes(ctx, w, h) {
  const bounds = plotBounds(w, h);
  const zeroY = yToPixel(0, w, h);
  const xTicks = [-100, -75, -50, -25, 0, 25, 50, 75, 100];
  const yTicks = [-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2];
  ctx.save();
  ctx.strokeStyle = "rgba(237,244,242,0.55)";
  ctx.fillStyle = "rgba(237,244,242,0.78)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(bounds.left, bounds.top);
  ctx.lineTo(bounds.left, bounds.bottom);
  ctx.moveTo(bounds.left, zeroY);
  ctx.lineTo(bounds.right, zeroY);
  ctx.stroke();

  ctx.font = "14px ui-sans-serif, system-ui, sans-serif";
  ctx.textAlign = "center";
  xTicks.forEach((tick) => {
    const px = xToPixel(tick, w, h);
    ctx.beginPath();
    ctx.moveTo(px, bounds.bottom);
    ctx.lineTo(px, bounds.bottom + 5);
    ctx.stroke();
    ctx.fillText(String(tick), px, bounds.bottom + 22);
  });

  ctx.textAlign = "right";
  yTicks.forEach((tick) => {
    const py = yToPixel(tick, w, h);
    ctx.beginPath();
    ctx.moveTo(bounds.left - 5, py);
    ctx.lineTo(bounds.left, py);
    ctx.stroke();
    ctx.fillText(formatTick(tick), bounds.left - 9, py + 4);
  });

  ctx.font = "16px ui-sans-serif, system-ui, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText("x", bounds.right - 10, zeroY - 10);
  ctx.translate(18, h * 0.5);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("psi(x)", 0, 0);
  ctx.restore();
}

function formatTick(valueAtTick) {
  if (Number.isInteger(valueAtTick)) return String(valueAtTick);
  const fixed = valueAtTick.toFixed(1);
  if (valueAtTick > 0 && valueAtTick < 1) return fixed.slice(1);
  if (valueAtTick < 0 && valueAtTick > -1) return `-${fixed.slice(2)}`;
  return fixed;
}

function drawPotentialOverlay(ctx, w, h, maxV) {
  ctx.save();
  ctx.strokeStyle = "#f0c36a";
  ctx.lineWidth = 2;
  ctx.setLineDash([7, 5]);
  ctx.beginPath();
  for (let i = 0; i < N; i += 1) {
    const px = xToPixel(x[i], w, h);
    const normalizedV = state.potential[i] / maxV;
    const py = yToPixel(normalizedV, w, h);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#f0c36a";
  ctx.font = "14px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("V(x)/Vmax", w - 132, 32);
  ctx.restore();
}

function potentialDisplayScale() {
  if (state.mode === "tunnel") return Number(controls.barrierHeight.max);
  return Math.max(...state.potential, 1e-12);
}

function drawWaveCurve(ctx, values, color, w, h, thicker = false) {
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = thicker ? 3 : 2;
  ctx.beginPath();
  values.forEach((valueAtX, i) => {
    const px = xToPixel(x[i], w, h);
    const py = yToPixel(valueAtX, w, h);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  });
  ctx.stroke();
  ctx.restore();
}

function drawWaveReadout(ctx, m, w) {
  ctx.save();
  ctx.fillStyle = "rgba(13,17,16,0.72)";
  ctx.fillRect(56, 22, 278, 100);
  ctx.strokeStyle = "rgba(255,255,255,0.12)";
  ctx.strokeRect(56, 22, 278, 100);
  ctx.fillStyle = "#edf4f2";
  ctx.font = "15px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText(`t = ${state.t.toFixed(3)}`, 70, 46);
  ctx.fillText(`Normalization = ${m.norm.toFixed(6)}`, 70, 68);
  ctx.fillText(`Energy = ${m.energy.toFixed(6)}`, 70, 90);
  ctx.fillText(`Width Delta x = ${m.width.toFixed(3)}`, 70, 112);
  ctx.fillStyle = "rgba(237,244,242,0.58)";
  ctx.fillText("time varying psi(x,t)", w - 188, 58);
  ctx.restore();
}

function drawWidthMarker(ctx, m, w, h) {
  const left = Math.max(-L, m.meanX - m.width);
  const right = Math.min(L, m.meanX + m.width);
  const y = yToPixel(0.08, w, h);
  const leftPx = xToPixel(left, w, h);
  const rightPx = xToPixel(right, w, h);
  ctx.save();
  ctx.strokeStyle = "rgba(127,181,255,0.72)";
  ctx.fillStyle = "rgba(127,181,255,0.9)";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(leftPx, y);
  ctx.lineTo(rightPx, y);
  ctx.stroke();
  ctx.beginPath();
  ctx.moveTo(leftPx, y - 6);
  ctx.lineTo(leftPx, y + 6);
  ctx.moveTo(rightPx, y - 6);
  ctx.lineTo(rightPx, y + 6);
  ctx.stroke();
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("Delta x", Math.min(rightPx + 8, w - 70), y + 4);
  ctx.restore();
}

function drawProbabilityDensity() {
  const canvas = els.probabilityCanvas;
  const ctx = canvas.getContext("2d");
  const w = canvas.width;
  const h = canvas.height;
  const density = Array.from({ length: N }, (_, i) => state.psiRe[i] ** 2 + state.psiIm[i] ** 2);
  const maxDensity = Math.max(...density, 1e-12);
  ctx.clearRect(0, 0, w, h);
  drawGrid(ctx, w, h);
  drawProbabilityAxes(ctx, w, h, maxDensity);

  ctx.save();
  ctx.strokeStyle = "#74d3ae";
  ctx.fillStyle = "rgba(116,211,174,0.18)";
  ctx.lineWidth = 2;
  ctx.beginPath();
  density.forEach((rho, i) => {
    const px = probabilityXToPixel(x[i], w, h);
    const py = probabilityYToPixel(rho, maxDensity, w, h);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  });
  ctx.stroke();
  ctx.lineTo(probabilityXToPixel(L, w, h), probabilityYToPixel(0, maxDensity, w, h));
  ctx.lineTo(probabilityXToPixel(-L, w, h), probabilityYToPixel(0, maxDensity, w, h));
  ctx.closePath();
  ctx.fill();
  ctx.fillStyle = "rgba(237,244,242,0.74)";
  ctx.font = "14px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("|psi(x)|^2", w - 104, 28);
  ctx.restore();
}

function probabilityBounds(w, h) {
  return {
    left: 54,
    right: w - 22,
    top: 20,
    bottom: h - 34,
  };
}

function probabilityXToPixel(xValue, w, h) {
  const bounds = probabilityBounds(w, h);
  return bounds.left + ((xValue + L) / (2 * L)) * (bounds.right - bounds.left);
}

function probabilityYToPixel(valueAtX, maxDensity, w, h) {
  const bounds = probabilityBounds(w, h);
  return bounds.bottom - (valueAtX / maxDensity) * (bounds.bottom - bounds.top);
}

function drawProbabilityAxes(ctx, w, h, maxDensity) {
  const bounds = probabilityBounds(w, h);
  const xTicks = [-100, -50, 0, 50, 100];
  ctx.save();
  ctx.strokeStyle = "rgba(237,244,242,0.5)";
  ctx.fillStyle = "rgba(237,244,242,0.72)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(bounds.left, bounds.top);
  ctx.lineTo(bounds.left, bounds.bottom);
  ctx.lineTo(bounds.right, bounds.bottom);
  ctx.stroke();
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.textAlign = "center";
  xTicks.forEach((tick) => {
    const px = probabilityXToPixel(tick, w, h);
    ctx.beginPath();
    ctx.moveTo(px, bounds.bottom);
    ctx.lineTo(px, bounds.bottom + 5);
    ctx.stroke();
    ctx.fillText(String(tick), px, bounds.bottom + 20);
  });
  ctx.textAlign = "right";
  ctx.fillText("0", bounds.left - 8, bounds.bottom + 4);
  ctx.fillText(maxDensity.toFixed(3), bounds.left - 8, bounds.top + 4);
  ctx.textAlign = "left";
  ctx.fillText("x", bounds.right - 10, bounds.bottom - 8);
  ctx.restore();
}

function drawEnergyPosition(m) {
  const canvas = els.energyPositionCanvas;
  const ctx = canvas.getContext("2d");
  const w = canvas.width;
  const h = canvas.height;
  const yMax = energyPlotMax(m.energy);
  ctx.clearRect(0, 0, w, h);
  drawGrid(ctx, w, h);
  drawEnergyAxes(ctx, w, h, yMax);
  drawEnergyPotential(ctx, w, h, yMax);
  drawEnergyLevel(ctx, w, h, yMax, m.energy);
  drawMeanPositionMarker(ctx, w, h, m.meanX);
}

function energyPlotMax(energyValue) {
  if (state.mode === "harmonic") return Math.max(Math.abs(energyValue) * 1.35, 1);
  if (state.mode === "tunnel") return Math.max(Number(controls.barrierHeight.max), Math.abs(energyValue) * 1.2, 1);
  return Math.max(Math.abs(energyValue) * 1.4, 1);
}

function energyBounds(w, h) {
  return {
    left: 54,
    right: w - 22,
    top: 20,
    bottom: h - 34,
  };
}

function energyXToPixel(xValue, w, h) {
  const bounds = energyBounds(w, h);
  return bounds.left + ((xValue + L) / (2 * L)) * (bounds.right - bounds.left);
}

function energyYToPixel(energyValue, yMax, w, h) {
  const bounds = energyBounds(w, h);
  const clamped = Math.max(0, Math.min(energyValue, yMax));
  return bounds.bottom - (clamped / yMax) * (bounds.bottom - bounds.top);
}

function drawEnergyAxes(ctx, w, h, yMax) {
  const bounds = energyBounds(w, h);
  const xTicks = [-100, -50, 0, 50, 100];
  ctx.save();
  ctx.strokeStyle = "rgba(237,244,242,0.5)";
  ctx.fillStyle = "rgba(237,244,242,0.72)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(bounds.left, bounds.top);
  ctx.lineTo(bounds.left, bounds.bottom);
  ctx.lineTo(bounds.right, bounds.bottom);
  ctx.stroke();
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.textAlign = "center";
  xTicks.forEach((tick) => {
    const px = energyXToPixel(tick, w, h);
    ctx.beginPath();
    ctx.moveTo(px, bounds.bottom);
    ctx.lineTo(px, bounds.bottom + 5);
    ctx.stroke();
    ctx.fillText(String(tick), px, bounds.bottom + 20);
  });
  ctx.textAlign = "right";
  ctx.fillText("0", bounds.left - 8, bounds.bottom + 4);
  ctx.fillText(yMax.toFixed(2), bounds.left - 8, bounds.top + 4);
  ctx.textAlign = "left";
  ctx.fillText("x", bounds.right - 10, bounds.bottom - 8);
  ctx.restore();
}

function drawEnergyPotential(ctx, w, h, yMax) {
  ctx.save();
  ctx.strokeStyle = "#f0c36a";
  ctx.lineWidth = 2;
  ctx.beginPath();
  let drawing = false;
  for (let i = 0; i < N; i += 1) {
    const px = energyXToPixel(x[i], w, h);
    const visiblePotential = Math.min(state.potential[i], yMax);
    const py = energyYToPixel(visiblePotential, yMax, w, h);
    if (!drawing) {
      ctx.moveTo(px, py);
      drawing = true;
    }
    else ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.fillStyle = "#f0c36a";
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("V(x)", w - 72, 27);
  ctx.restore();
}

function drawEnergyLevel(ctx, w, h, yMax, energyValue) {
  const bounds = energyBounds(w, h);
  const py = energyYToPixel(energyValue, yMax, w, h);
  ctx.save();
  ctx.strokeStyle = "#7fb5ff";
  ctx.fillStyle = "#7fb5ff";
  ctx.lineWidth = 2;
  ctx.setLineDash([7, 5]);
  ctx.beginPath();
  ctx.moveTo(bounds.left, py);
  ctx.lineTo(bounds.right, py);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText(`<E> = ${energyValue.toFixed(3)}`, bounds.left + 8, Math.max(py - 8, bounds.top + 14));
  ctx.restore();
}

function drawMeanPositionMarker(ctx, w, h, meanX) {
  const bounds = energyBounds(w, h);
  const px = energyXToPixel(meanX, w, h);
  ctx.save();
  ctx.strokeStyle = "rgba(116,211,174,0.75)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(px, bounds.top);
  ctx.lineTo(px, bounds.bottom);
  ctx.stroke();
  ctx.fillStyle = "rgba(116,211,174,0.9)";
  ctx.font = "13px ui-sans-serif, system-ui, sans-serif";
  ctx.fillText("<x>", Math.min(px + 5, bounds.right - 28), bounds.top + 16);
  ctx.restore();
}

function drawGrid(ctx, w, h) {
  ctx.fillStyle = "#0d1110";
  ctx.fillRect(0, 0, w, h);
  ctx.strokeStyle = "rgba(255,255,255,0.08)";
  ctx.lineWidth = 1;
  for (let i = 1; i < 5; i += 1) {
    const y = (i / 5) * h;
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(w, y);
    ctx.stroke();
  }
}

function updateReadouts(measurement) {
  const m = measurement || measure();
  const normError = Math.abs(1 - m.norm);
  const highK = Math.abs(value("k0")) + 3 / value("sigma") > Math.max(...k.map(Math.abs)) * 0.72;
  const warning = normError > 0.02 || highK || value("dt") > 0.032;
  els.timeReadout.textContent = state.t.toFixed(3);
  els.normReadout.textContent = m.norm.toFixed(4);
  els.energyReadout.textContent = m.energy.toFixed(3);
  els.meanXReadout.textContent = m.meanX.toFixed(3);
  els.meanPReadout.textContent = m.meanP.toFixed(3);
  els.widthReadout.textContent = m.width.toFixed(3);
  els.stripTime.textContent = state.t.toFixed(3);
  els.stripNorm.textContent = m.norm.toFixed(4);
  els.stripEnergy.textContent = m.energy.toFixed(3);
  els.stripMeanX.textContent = m.meanX.toFixed(3);
  els.stripMeanP.textContent = m.meanP.toFixed(3);
  els.stripWidth.textContent = m.width.toFixed(3);
  els.alertReadout.textContent = warning ? (highK ? "Aliasing risk" : "Check dt/norm") : "Stable";
  els.alertReadout.classList.toggle("is-warning", warning);
}

function render() {
  resizeCanvases();
  const m = measure();
  drawWave();
  drawProbabilityDensity();
  drawEnergyPosition(m);
  updateReadouts(m);
}

function resizeCanvases() {
  [els.waveCanvas, els.probabilityCanvas, els.energyPositionCanvas].forEach((canvas) => {
    const rect = canvas.getBoundingClientRect();
    const nextWidth = Math.max(320, Math.round(rect.width));
    const nextHeight = Math.max(120, Math.round(rect.height));
    if (canvas.width !== nextWidth) canvas.width = nextWidth;
    if (canvas.height !== nextHeight) canvas.height = nextHeight;
  });
}

function reset() {
  updatePotential();
  initializeWavepacket();
  render();
}

function frame() {
  if (state.running) {
    const steps = value("speed");
    for (let i = 0; i < steps; i += 1) step();
    render();
  }
  requestAnimationFrame(frame);
}

document.querySelectorAll(".mode-tab").forEach((tab) => {
  tab.addEventListener("click", () => {
    document.querySelectorAll(".mode-tab").forEach((t) => t.classList.remove("active"));
    tab.classList.add("active");
    state.mode = tab.dataset.mode;
    setDefaultsForMode(state.mode);
    reset();
  });
});

Object.entries(controls).forEach(([id, control]) => {
  control.addEventListener("input", () => {
    updateOutputs();
    if (["x0", "k0", "sigma"].includes(id)) {
      reset();
      return;
    }
    if (["omega", "barrierHeight", "barrierWidth"].includes(id)) updatePotential();
    render();
  });
});

els.playBtn.addEventListener("click", () => {
  state.running = !state.running;
  els.playBtn.textContent = state.running ? "Pause" : "Play";
});

els.resetBtn.addEventListener("click", reset);
window.addEventListener("resize", render);

updateOutputs();
reset();
requestAnimationFrame(frame);

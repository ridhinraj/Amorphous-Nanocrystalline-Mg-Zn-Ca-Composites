const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.defineLayout({ name: "FIG1", width: 13.333, height: 6.10 });
pres.layout = "FIG1";

// ---- palette ----
const SLATE = "3B5998";
const CHARCOAL = "2B2D42";
const TEAL = "006D77";
const ORANGE = "D95F02";
const GREEN = "1B9E77";
const MAGENTA = "C51B8A";
const GOLD = "E6AB02";
const GOLD_TEXT = "CC9600";
const BODY = "45453F";
const REF = "8A8378";
const FOOTER_BG = "F4F3EE";
const FOOTER_BORDER = "E8E6DD";
const CARD_BORDER = "E3E1D9";
const TITLE_COL = "22242B";
const WHITE = "FFFFFF";

const slide = pres.addSlide();
slide.background = { color: WHITE };

function card(opts) {
  const { x, y, w, h, accent } = opts;
  slide.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.07,
    fill: { color: WHITE },
    line: { color: CARD_BORDER, width: 0.75 },
  });
  // top accent bar
  slide.addShape(pres.ShapeType.rect, {
    x, y, w, h: 0.045,
    fill: { color: accent }, line: { type: "none" },
  });
}

function iconCircle(opts) {
  const { x, y, letter, color } = opts;
  slide.addShape(pres.ShapeType.ellipse, {
    x, y, w: 0.34, h: 0.34,
    fill: { color }, line: { type: "none" },
  });
  slide.addText(letter, {
    x, y, w: 0.34, h: 0.34, margin: 0,
    align: "center", valign: "middle",
    fontFace: "Arial", fontSize: 11, bold: true, color: WHITE,
    isTextBox: true,
  });
}

function cardHead(opts) {
  const { x, y, w, letter, accent, title, ref } = opts;
  iconCircle({ x: x + 0.14, y: y + 0.12, letter, color: accent });
  slide.addText(title, {
    x: x + 0.56, y: y + 0.10, w: w - 0.7, h: 0.24, margin: 0,
    fontFace: "Arial", fontSize: 12.5, bold: true, color: TITLE_COL, valign: "bottom", isTextBox: true,
  });
  slide.addText(ref, {
    x: x + 0.56, y: y + 0.34, w: w - 0.7, h: 0.16, margin: 0,
    fontFace: "Arial", fontSize: 8.5, color: REF, valign: "top", isTextBox: true,
  });
}

function bulletsBox(opts) {
  const { x, y, w, h, items, fontSize } = opts;
  const fs = fontSize || 10.5;
  const finalRuns = [];
  items.forEach((it) => {
    const runs = Array.isArray(it) ? it : [{ text: it, options: {} }];
    runs.forEach((r, j) => {
      const text = j === 0 ? "–  " + r.text : r.text;
      finalRuns.push({
        text,
        options: Object.assign(
          { fontFace: "Arial", fontSize: fs, color: BODY, breakLine: j === runs.length - 1 },
          r.options || {}
        ),
      });
    });
  });
  slide.addText(finalRuns, {
    x, y, w, h, margin: 0, valign: "top",
    lineSpacingMultiple: 1.12, paraSpaceAfter: 3, isTextBox: true,
  });
}

function footerStrip(opts) {
  const { x, y, w, h, label, text } = opts;
  slide.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.05,
    fill: { color: FOOTER_BG }, line: { color: FOOTER_BORDER, width: 0.75 },
  });
  slide.addText(
    [
      { text: label + "  ", options: { bold: true, color: CHARCOAL, fontSize: 9.5 } },
      { text: text, options: { color: "5C5C57", fontSize: 9.5 } },
    ],
    { x: x + 0.14, y, w: w - 0.28, h, margin: 0, valign: "middle", fontFace: "Arial", isTextBox: true }
  );
}

function arrow(x1, y1, x2, y2, opts) {
  const o = Object.assign({ color: CHARCOAL, width: 1.5, dash: "solid" }, opts || {});
  slide.addShape(pres.ShapeType.line, {
    x: Math.min(x1, x2), y: Math.min(y1, y2),
    w: Math.abs(x2 - x1) || 0.01, h: Math.abs(y2 - y1) || 0.01,
    flipV: y2 < y1 ? false : false,
    line: {
      color: o.color, width: o.width, dashType: o.dash,
      endArrowType: "triangle",
    },
  });
}


const SUP = (t) => ({ text: t, options: { superscript: true } });
const SUB = (t) => ({ text: t, options: { subscript: true } });
const T = (t) => ({ text: t, options: {} });
const COLS = [0.427, 3.613, 6.800, 9.987], CW = 2.920;

// ================= Row (a) =================
slide.addText("(a) Model construction", { x: 0.427, y: 0.22, w: 6, h: 0.3, margin: 0,
  fontFace: "Arial", fontSize: 15, bold: true, color: SLATE, isTextBox: true });
const AY = 0.58, AH = 1.80;
const A = [
  ["P", "Interatomic potential", "Sec. 2.1, 3.1", [
     "2NN MEAM Mg–Zn–Ca + ZBL below 2 Å, in LAMMPS",
     "Checked: lattice, elastic constants, basal stacking faults, mixing enthalpy"]],
  ["M", "Amorphous matrix", "Sec. 2.2, 2.6, 3.2", [
     "389,344 atoms, ~20 nm cubic cell, 10 compositions",
     [T("Melt at 1500 K, quench at 10"), SUP("12"), T(" K/s, 600 K cycle")],
     "RDF, Tg and fictive temperature"]],
  ["C", "Composite construction", "Sec. 2.3", [
     "Pure-glass controls (no seed)",
     "One HCP α-Mg seed, R = 15–35 Å",
     "12–20 Å seeds at 5–19% volume fraction, one orientation",
     "Four packings per matrix at 10%"]],
  ["T", "Relaxation and tension", "Sec. 2.4, 2.5", [
     "Minimization, then 50 ps NPT at 300 K",
     [T("Uniaxial tension to 15% at 10"), SUP("9"), T(" s"), SUP("−1"), T(", lateral stress zero")],
     [T("σ"), SUB("max"), T(" (peak strength), E, σ"), SUB("15"), T(", softening")]]],
];
A.forEach(([L, title, ref, items], i) => {
  card({ x: COLS[i], y: AY, w: CW, h: AH, accent: SLATE });
  cardHead({ x: COLS[i], y: AY, w: CW, letter: L, accent: SLATE, title, ref });
  bulletsBox({ x: COLS[i] + 0.14, y: AY + 0.56, w: CW - 0.26, h: AH - 0.62, items, fontSize: 10 });
});

arrow(6.667, AY + AH, 6.667, 2.68, { color: CHARCOAL, width: 1.75 });

// ================= Row (b) =================
const BY = 2.68;
card({ x: 0.427, y: BY, w: 12.480, h: 1.10, accent: CHARCOAL });
slide.addText([
    { text: "(b) 48 tensile simulations", options: { bold: true, fontSize: 15, color: TITLE_COL } },
    { text: "   Sec. 3.3–3.7 · first campaign 25 runs, second campaign 23 runs", options: { fontSize: 9.5, color: REF } } ],
  { x: 0.56, y: BY + 0.07, w: 12.2, h: 0.32, margin: 0, valign: "middle", fontFace: "Arial", isTextBox: true });
function pill(x, color, title, detail) {
  const y = BY + 0.48, w = 2.93, h = 0.50;
  slide.addShape(pres.ShapeType.roundRect, { x, y, w, h, rectRadius: 0.06, fill: { color, transparency: 90 }, line: { type: "none" } });
  slide.addShape(pres.ShapeType.rect, { x, y, w: 0.045, h, fill: { color }, line: { type: "none" } });
  slide.addText([
      { text: title, options: { bold: true, fontSize: 11, color, breakLine: true } },
      { text: detail, options: { fontSize: 9.5, color: "5C5C57" } } ],
    { x: x + 0.16, y, w: w - 0.26, h, margin: 0, valign: "middle", fontFace: "Arial", isTextBox: true });
}
pill(0.577, ORANGE, "Composition", "10 alloys, 3 pure glasses, 2 HCP crystals");
pill(3.657, GREEN, "Volume fraction", "5–19% in three matrices");
pill(6.737, MAGENTA, "Seed radius", "12, 15 and 20 Å at 10%");
pill(9.817, GOLD_TEXT, "Rate and temperature", "10⁸–2×10⁹ s⁻¹; 250–400 K");

arrow(6.667, BY + 1.10, 6.667, 4.25, { color: CHARCOAL, width: 1.75 });
arrow(11.447, 4.25, 11.447, BY + 1.10, { color: REF, width: 1.5, dash: "dash" });
slide.addText("ranks follow-up runs", { x: 11.55, y: 3.86, w: 1.4, h: 0.3, margin: 0,
  fontFace: "Arial", fontSize: 9, italic: true, color: REF, isTextBox: true });

// ================= Row (c) =================
slide.addText("(c) Analysis", { x: 0.427, y: 3.90, w: 6, h: 0.3, margin: 0,
  fontFace: "Arial", fontSize: 15, bold: true, color: TEAL, isTextBox: true });
const CY = 4.25, CH = 1.62;
const C = [
  ["Σ", "Bulk mechanics", "Sec. 3.3, 3.4, 3.7", [
     [T("σ"), SUB("max"), T(", E, flow stress, softening")],
     "Against composition, volume fraction, seed radius, rate and temperature"]],
  ["D", "Plasticity and interphase", "Sec. 2.7–2.9, 3.5", [
     "Von Mises stress test for a shear band",
     [T("Incremental D"), SUP("2"), T("min clusters (DBSCAN)")],
     "3 Å interphase shell vs. bulk, 236 seeds"]],
  ["V", "Ca and Zn environments", "Sec. 2.10, 2.11, 3.6", [
     "Voronoi volume-normalized stress",
     [T("q"), SUB("4"), T(", q"), SUB("6"), T(" and Warren–Cowley α")],
     "Gradient boosting with SHAP"]],
  ["GP", "Surrogate model", "Sec. 2.12, 3.8", [
     "Matérn kernel, trained on peak strength",
     "Uncertainty ranks follow-up runs",
     "Checked on later, held-out runs"]],
];
C.forEach(([L, title, ref, items], i) => {
  card({ x: COLS[i], y: CY, w: CW, h: CH, accent: TEAL });
  cardHead({ x: COLS[i], y: CY, w: CW, letter: L, accent: TEAL, title, ref });
  bulletsBox({ x: COLS[i] + 0.14, y: CY + 0.56, w: CW - 0.26, h: CH - 0.62, items, fontSize: 10 });
});

pres.writeFile({ fileName: "/tmp/w/fig15/f1/Figure_01_schematic.pptx" }).then(() => console.log("done"));

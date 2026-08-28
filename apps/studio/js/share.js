// A run, in a URL fragment.
//
// The point of this project is that a run is a pure function of its scenario:
// the same description gives bit-identical output on x86-64, aarch64 and wasm.
// That claim is only useful if the description travels. A link is the cheapest
// vehicle there is, and it is what lets an instructor hand out an exact
// scenario rather than a data file, and what lets a bug report say "here" and
// mean it.
//
// # The fragment, not the query string
//
// Everything after `#` stays in the browser and is never sent to the server.
// GitHub Pages and a Hugging Face Static Space serve one static page, so a
// query string would be an unused round trip carrying a scenario into somebody
// else's logs. The fragment also changes without a navigation, so the address
// bar can track the panel as it is edited.
//
// # Defaults are omitted
//
// A link for the baseline is `#`. Only fields that differ from the defaults are
// written, so the common link is short and a reader can see at a glance what
// was changed. That also means a future field with a sensible default does not
// invalidate the links already in circulation: an old link simply does not
// mention it and gets the default.
//
// # Nothing here trusts its input
//
// A fragment is user input from an untrusted link. Every field is range checked
// against what the bindings will accept, and anything unparseable falls back to
// the default rather than throwing, because a mistyped link should open the
// baseline plant and not a blank page. The bindings validate again on their own
// side, which is where a rejection would actually be produced.

/** Short keys, so a link with a fault set and a channel list stays readable. */
const KEYS = {
  seed: "s",
  hours: "h",
  stepHours: "dt",
  sampleEvery: "n",
  integrator: "i",
  controlled: "c",
  driverForcesIdv12: "d",
  tripEndsTheRun: "t",
  faults: "f",
  channels: "p",
  chunkSamples: "k",
  speedMultiple: "x",
};

const INTEGRATORS = new Set(["euler", "rk4", "dopri5", "dormand-prince"]);

function finiteOr(text, fallback, { min = -Infinity, max = Infinity } = {}) {
  const value = Number(text);
  if (!Number.isFinite(value) || value < min || value > max) return fallback;
  return value;
}

function boolOr(text, fallback) {
  if (text === "1") return true;
  if (text === "0") return false;
  return fallback;
}

/** A comma-separated list of integers, filtered to a range and de-duplicated. */
function intListOr(text, fallback, min, max) {
  if (text === undefined) return fallback;
  if (text === "") return [];
  const seen = new Set();
  for (const part of text.split(",")) {
    const value = Number(part);
    if (!Number.isInteger(value) || value < min || value > max) continue;
    seen.add(value);
  }
  return [...seen].sort((a, b) => a - b);
}

/**
 * Render the parts of `state` that differ from `defaults` as a fragment.
 *
 * The leading `#` is not included, so the caller decides whether an empty
 * result means `#` or an untouched address bar.
 *
 * @param {object} state
 * @param {object} defaults
 * @returns {string}
 */
export function encodeState(state, defaults) {
  const parts = [];
  const put = (key, value) => parts.push(`${key}=${value}`);

  for (const [field, key] of Object.entries(KEYS)) {
    const value = state[field];
    const fallback = defaults[field];
    if (value === undefined) continue;

    if (Array.isArray(value)) {
      const a = [...value].sort((x, y) => x - y).join(",");
      const b = [...(fallback ?? [])].sort((x, y) => x - y).join(",");
      if (a !== b) put(key, a);
      continue;
    }
    if (typeof value === "boolean") {
      if (value !== fallback) put(key, value ? "1" : "0");
      continue;
    }
    if (typeof value === "number") {
      // The seed is 4651207995, an integer that must survive the round trip
      // exactly; `toString` on a double gives the shortest text that reparses
      // to the same double, which is precisely the guarantee needed.
      if (value !== fallback) put(key, String(value));
      continue;
    }
    if (value !== fallback) put(key, encodeURIComponent(String(value)));
  }

  return parts.join("&");
}

/**
 * Parse a fragment back into a state, falling back to `defaults` field by
 * field. Never throws.
 *
 * @param {string} fragment with or without a leading `#`
 * @param {object} defaults
 * @param {number} channelCount how many columns a packed row has
 * @returns {object}
 */
export function decodeState(fragment, defaults, channelCount = 54) {
  const state = { ...defaults, faults: [...defaults.faults], channels: [...defaults.channels] };
  const text = fragment.startsWith("#") ? fragment.slice(1) : fragment;
  if (text === "") return state;

  const raw = new Map();
  for (const pair of text.split("&")) {
    if (pair === "") continue;
    const eq = pair.indexOf("=");
    if (eq < 0) continue;
    raw.set(pair.slice(0, eq), decodeURIComponent(pair.slice(eq + 1)));
  }
  const get = (field) => raw.get(KEYS[field]);

  // The seed must be finite and positive; `teprob.f:1187` compiles in
  // 4651207995 and the generator is undefined for anything else.
  state.seed = finiteOr(get("seed"), defaults.seed, { min: Number.MIN_VALUE });
  state.hours = finiteOr(get("hours"), defaults.hours, { min: Number.MIN_VALUE });
  state.stepHours = finiteOr(get("stepHours"), defaults.stepHours, {
    min: Number.MIN_VALUE,
  });

  const cadence = finiteOr(get("sampleEvery"), defaults.sampleEvery, { min: 1 });
  state.sampleEvery = Math.max(1, Math.round(cadence));

  const chunk = finiteOr(get("chunkSamples"), defaults.chunkSamples, { min: 1 });
  state.chunkSamples = Math.max(1, Math.round(chunk));

  // Zero is a legal speed and means "as fast as this machine goes".
  state.speedMultiple = finiteOr(get("speedMultiple"), defaults.speedMultiple, {
    min: 0,
  });

  const integrator = get("integrator");
  state.integrator = INTEGRATORS.has(integrator) ? integrator : defaults.integrator;

  state.controlled = boolOr(get("controlled"), defaults.controlled);
  state.driverForcesIdv12 = boolOr(
    get("driverForcesIdv12"),
    defaults.driverForcesIdv12,
  );
  state.tripEndsTheRun = boolOr(get("tripEndsTheRun"), defaults.tripEndsTheRun);

  // Twenty, not the twenty-one of the later literature: `teprob.f:340` loops
  // `DO 500 I=1,20`.
  state.faults = intListOr(get("faults"), defaults.faults, 1, 20);
  // Column 0 is the time and is not plottable.
  state.channels = intListOr(get("channels"), defaults.channels, 1, channelCount - 1);

  return state;
}

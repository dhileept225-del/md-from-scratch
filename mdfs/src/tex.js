import katex from 'katex';
import 'katex/dist/katex.min.css';

export const tex = (s, display = false) =>
  katex.renderToString(s, { displayMode: display, throwOnError: false });

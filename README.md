# Quantum Wavepacket Lab

A browser-based 1D quantum mechanics simulator for exploring:

- Quantum particle motion in a harmonic trap
- Quantum tunneling through a finite barrier
- Free Gaussian wavepacket spreading

The UI uses a split-step Fourier propagation loop with visualizations for time-varying `Re psi(x)`, `Im psi(x)`, `|psi(x)|`, probability density `|psi(x)|^2`, normalized `V(x)/Vmax`, energy versus position, normalization, and key diagnostics. The main plot uses `x = [-100, 100]` and `psi(x) = [-0.4, 1.2]`. Diagnostics include energy, `<x>`, and `<p>`. The harmonic oscillator defaults are tuned for visible oscillatory motion, and the layout is compact so the main plots and controls can fit into a single screenshot on typical desktop screens.

Open `index.html` in a browser to run the simulator on desktop.

For Android sharing/opening through WhatsApp or a file picker, use `Quantum_Wavepacket_Lab_Android.html`. It is a single self-contained HTML file with CSS and JavaScript embedded, so Android does not need to load separate `styles.css` or `app.js` files.

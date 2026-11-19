# Monica Dessole — personal website

A minimal, responsive academic website, published with GitHub Pages.
The homepage uses plain HTML, CSS, and a small script for publication filtering;
it also works without JavaScript or externally hosted fonts and UI libraries.

## Local preview

From this directory, run:

```sh
python3 -m http.server 8000
```

Open <http://localhost:8000>. Stop the server with Ctrl+C.
The existing Jekyll configuration is retained for GitHub Pages.

## Updating the site

- Edit biography, links, software, teaching, and publications in `index.html`.
- Publication entries use `data-type="journal"`, `"conference"`, `"preprint"`, or
  `"thesis"`. The filtering controls are in the same file.
- Adjust colors, spacing, and responsive layouts in `assets/style.css`.
- Replace the portrait in `images/IMG_9995.JPG` and the CV in `pdfs/cv_dessole.pdf`
  when needed.

Publication links open in the same tab; the CV and thesis link directly to PDFs.
All publications remain visible when JavaScript is disabled and when printing.

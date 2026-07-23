# Communications

Abstracts, posters, slides, and other material written **about** `synpmx` for an
outside audience. This is the counterpart to `references/`, which holds material
written by other people that the package draws on.

Not shipped with the package: `.Rbuildignore` excludes this folder, so nothing
here reaches a user who installs `synpmx`.

## What goes here

Keep the source text and the rendered artifact side by side, named by year and
venue:

```
2026-acop-abstract.md      the text, in Markdown
2026-acop-poster.pdf       the rendered artifact
```

The Markdown source matters as much as the PDF. Abstracts get revised, rejected,
retargeted at another venue, and mined for README and vignette prose. A
diffable source makes all of that tractable; a PDF makes none of it.

## Before anything here becomes public

Two checks, and neither is a formality given what this package is for.

- **No confidential data.** Figures, tables, and worked examples must come from
  public datasets or package fixtures. A plot built from an internal study is
  source-derived output and belongs in `scripts_private/`, whatever it looks
  like once rendered.
- **Internal review first.** Clear the abstract through whatever review your
  organization requires before submission, not after.

## Publishing to the website

`pkgdown/assets/` is copied verbatim to the site root, so a file placed there is
served at a stable public URL — useful as a QR-code target on a poster. Copy a
poster there only once both checks above have passed. Putting it in this folder
publishes nothing; putting it in `pkgdown/assets/` publishes it to anyone.

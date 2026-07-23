# Capture a schema asserted to be public

This helper records names and practical classes only. Calling it is an
assertion that those metadata and factor/category levels are public; do
not derive them from confidential data unless their release has been
approved.

## Usage

``` r
pmx_schema(data, exclude = NULL)
```

## Arguments

- data:

  A public schema template.

- exclude:

  Columns omitted from generated output.

## Value

A `pmx_schema` object.

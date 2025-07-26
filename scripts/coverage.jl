using LocalCoverage

coverage = generate_coverage("FLORIDyn"; run_test = true)
html_coverage(coverage; open = true, dir = tempdir())
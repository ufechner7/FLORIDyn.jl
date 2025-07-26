using LocalCoverage

coverage = generate_coverage("FLORIDyn"; run_test = true)
html_coverage(coverage; open = true, dir = tempdir())

@info "On Ubuntu 2024 you need to execute:"
@info "cd /tmp"
@info "python3 -m http.server 8000"
@info "And then open http://localhost:8000 with your browser."
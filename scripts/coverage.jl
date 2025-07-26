using LocalCoverage
import LiveServer as LS

coverage = generate_coverage("FLORIDyn"; run_test = true)
html_coverage(coverage; open = false, dir = tempdir())
LS.serve(launch_browser=true, dir="/tmp")

@info "On Ubuntu 2024 you need to execute:"
@info "cd /tmp"
@info "python3 -m http.server 8000"
@info "And then open http://localhost:8000 with your browser."
.PHONY: run monte-carlo test report clean all

run:
	python src/simulate_counts.py
	python src/threshold_analysis.py

monte-carlo:
	python src/threshold_analysis.py --n-trials 10000

test:
	python -m unittest discover -s tests

report:
	cd report && pdflatex -interaction=nonstopmode detector_threshold_demo_report.tex
	cd report && pdflatex -interaction=nonstopmode detector_threshold_demo_report.tex

all: run monte-carlo test report

clean:
	rm -rf src/__pycache__ tests/__pycache__
	rm -f report/*.aux report/*.log report/*.out report/*.toc report/*.lof report/*.lot report/*.fls report/*.fdb_latexmk report/*.synctex.gz

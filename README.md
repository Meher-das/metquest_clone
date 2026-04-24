# MetQuest - For identifying efficient pathways to produce metabolites

Setup
=====
Make sure pip is installed on your system by using the command ```sudo apt install python3-pip```.

Create a virtual environment using the command ```python3 -m venv .venv```.

Source the virtual environment and run pip install to complete environment setup
```source .venv/bin/activate```,
```pip install -r requirements.txt```

To run the results from the case studies in papers 1 and 2 run the script,
```python3 test/original_validation.py```

Cases Studies
=============
1. The "1,4-Butanediol" Benchmarking (E. coli) | Paper: Yim et al. (2011), Nature Chemical Biology.

2. The "Artemisinic Acid" Precursor (S. cerevisiae) | Paper: Ro et al. (2006), Nature / Paddon et al. (2013), Nature.

Citation 
========
Ravikrishnan, A., Nasre, M. & Raman, K. Enumerating all possible biosynthetic pathways in metabolic networks. *Sci. Rep.* **8**:9932 (2018).

Barra, A. L. C., Dantas, L. O. C., Morão, L. G., Gutierrez, R. F., Polikarpov, I., Wrenger, C. & Nascimento, A. S. Essential metabolic routes as a way to ESKAPE from antibiotic resistance. Front. Public Health 8, 26 (2020).
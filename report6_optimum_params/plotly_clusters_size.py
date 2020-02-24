import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from itertools import groupby
import os

all_references=['extractedGTF', 'longIso_extractedGTF', 'nonoverlap_longIso_extractedGTF', 'noPseudo_nonoverlap_longIso_extractedGTF']

def fasta_iter(fasta_name):
    fh = open(fasta_name, 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq



REF1="extractedGTF_gencode.v33.transcripts.fa"
REF2="longIso_extractedGTF_gencode.v33.transcripts.fa"
REF3="nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa"
REF4="noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa"
REFS = [REF1, REF2, REF3, REF4]


histo = dict()
sims = [91, 93, 95, 97, 99]

for sim in sims:
    histo[f"{sim}%"] = list()


for REF_FASTA in REFS:
    for sim in sims:
        file = f"clusters_0.{sim}_cDBG_{REF_FASTA}.clstr"
        large_clusters = 0
        for header, seq in fasta_iter("cdhit_clusters/" + file):
            size = seq.count('nt')
            if size > 1:
                large_clusters += 1

        title = REF_FASTA.split('.')[0].split('_')[:-1]
        title = "_".join(title)
        histo[f"{sim}%"].append(large_clusters)

data = list()

for k,v in histo.items():
    data.append(go.Bar(name=f'{k}', x=all_references, y=v))


fig = go.Figure(data)


# Change the bar mode
fig.update_layout(barmode='group',)
plot(fig, filename = 'plotly_cdhit.html', auto_open=False)

# Change the bar mode
fig.update_layout(barmode='group', yaxis_type="log")
plot(fig, filename = 'plotly_cdhit_log.html', auto_open=False)
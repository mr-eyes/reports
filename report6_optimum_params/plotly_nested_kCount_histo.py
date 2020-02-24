import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import os

all_references=['extractedGTF', 'longIso_extractedGTF', 'nonoverlap_longIso_extractedGTF', 'noPseudo_nonoverlap_longIso']


kmers = [25, 31,41,51,75,81, 91]
kmers = list(map(str, kmers))

r_data = dict()
for k in kmers:
    r_data[k] = list()

refs = [
"extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv",
"noPseudo_nonoverlap_longIso_extractedGTF_gencode.v33.transcripts.fa_kmers_histogram.tsv"
]

for ref in refs:
    with open(ref , 'r') as tsv_reader:
        for line in tsv_reader:
            line = line.strip().split()
            kSize = line[0]
            uniq = int(line[1])
            title = os.path.basename(ref).split('_')[0] # .split('_')[:-1]
            r_data[kSize].append(uniq)

data = list()

for k,v in r_data.items():
    data.append(go.Bar(name=f'k{k}', x=all_references, y=v))


fig = go.Figure(data)


# Change the bar mode
fig.update_layout(barmode='group')
plot(fig, filename = 'plotly_histo.html', auto_open=False)

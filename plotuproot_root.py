import sys
import uproot # nuroot its little lib for root 
import mplhep as hep #lib for ploting 
import matplotlib.pyplot as plt #

def process (filename="AliXiCascade.root"): # an
    try:
        with uproot open (filename) as a file:     
            print("<b>Keys in file:</b>")
            print(file.classnames())


            container=file["MyTask/MyOutputContainer;1"]
            for i, obj in enumerate(container):
            print(f"Item {i}: {obj}")

            histogram=[
                "fHistMult", 
                "fHistPt", 
                "fHistdEdx", 
                "fHistInvMassLambdap", 
                "HistInvMassLambda", 
                "cos0", 
                "cos1", 
                "cos10", 
                "cos50", 
                "fHistCosThetalambdaproton"
            ]
            found_histograms = {}

            for obj in container:
                if hasattr(obj, "name") and obj.name in histogram
                print(f"<b>Found hsito:</b> {obj.name}")

            return histogram

except filenotfound:
    print(f"<b>Error:</b> File {filename} not found")
    return None
except Exception as e:
    print(f"</b> line:</b> {e}")
    return None

except nohistogram:
    if not histogram:
    print(f"</b> ho histo found")
    return

plt.style.use(hep.style.ALICE)  # or other styles like CMS, ALICE

hep.histplot(hist_mu)
hist_mult = hist_mu.to_hist()
hist_mult.plot()

# Access the histogram using full path
plt.hist(hist_mult, bins='auto', edgecolor='black')
plt.figure(figsize=(8, 6))
plt.xlabel("X", fontsize=14)
plt.ylabel("Y", fontsize=14)
plt.title("Multiplicity", fontsize=16)
plt.savefig("plot.png", dpi=300, bbox_inches='tight')


hist_mult.plot()
plt.show(block=True)   

import sys
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import m_e, c


def process_root_file(filename="coswindows.root"):
    try:
        
        with uproot.open(filename) as file:
           
            print("<b>Keys in file:</b>")
            print(file.classnames())
            container = file["MyTask/MyOutputContainer;1"]

            histograms_of_interest = [

                
                "fHistInvMassLambdap", 
                "fHistInvMassLambda", 
                "fHistPt",  
                "cos0", 
                "cos1", 
                "cos10", 
                "cos50",
                "fHistMult",
                "fHistEta",
                "fHistdEdx",
                "fHistptvsalpha",
                "fHistCosThetalambdaproton",
                "fHistPtP"

            ]

          
            found_histograms = {}

            
            for obj in container:
                if hasattr(obj, "name") and obj.name in histograms_of_interest:
                    found_histograms[obj.name] = obj
                    print(f"<b>Found histogram:</b> {obj.name}")

            return found_histograms

    except FileNotFoundError:
        print(f"<b>Error:</b> File {filename} not found.")
        return None
    except Exception as e:
        print(f"<b>An error occurred:</b> {e}")
        return None
        

def plot_histograms(histograms):
    if not histograms:
        print("<b>No histograms to plot</b>")
        return


    plt.style.use(hep.style.ALICE)
    
    histoname = {
        
        "cos0": "cos[-1;-0.5]",
        "cos1": "cos[-0.5,0]",
        "cos10": "cos[0;0.5]",
        "cos50": "cos[0.5,1]",
      
    }

    colors = ['blue', 'cyan', 'purple', 'magenta']
   
    fig = plt.figure(figsize=(12, 8)) 
    gs = fig.add_gridspec(1, 1)       
    ax1 = fig.add_subplot(gs[0, 0])

    for (name, hist), color in zip([(k, v) for k, v in histograms.items() if k.startswith('cos')], colors):
     
        try:
            hist_converted = hist.to_hist()
            hep.histplot(hist_converted, 
                     ax=ax1, 
                     label=histoname.get(name, name), 
                     color=color)
        except Exception as e:
            print(f"<b>Error plotting {name}:</b> {e}")
       

    ax1.set_title("InvMassLambdaP", fontsize=16)  
    ax1.set_xlabel("Gev/c2", fontsize=14)       
    ax1.set_ylabel("Entries", fontsize=14)     
    ax1.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout() 
    ax1.legend(title="CosTheta", loc='best')
    plt.savefig("cosine_windows_combined.png", dpi=300, bbox_inches='tight')   
    plt.show()
  
    for name, hist in histograms.items():
        if not name.startswith('cos'):
            try:
                
                h = hist.to_hist()
                
                if h.ndim == 2:
                    
                    fig = plt.figure(figsize=(10, 8))
                    ax = fig.add_subplot()
                    hep.hist2dplot(h, ax=ax)
                    
                    if name == "fHistdEdx":

                        
                        def bethe_bloch(p, m):
                            beta = p / np.sqrt(p**2 + m**2)
                            gamma = 1 / np.sqrt(1 - beta**2)
                            return (1/beta**2) * (np.log(2*0.511*beta**2*gamma**2 / (23.2e-6)) - beta**2)

                        
                        ax.set_xscale('log')
                        p = np.logspace(np.log10(0.1), np.log10(10), 100)

                        for m, label, color in [(0.139, 'Pion', 'red'), (0.494, 'Kaon', 'blue'), (0.938, 'Proton', 'green')]:
                            ax.plot(p, bethe_bloch(p, m), color=color, label=label, linewidth=2)
                        
                    ax.legend()
                    ax.set_title(name)
                    plt.show()

                
                else:
                   
                    fig = plt.figure(figsize=(10, 7))
                    ax = fig.add_subplot()
                    hep.histplot(h, ax=ax, color='black')
                    ax.set_title(name)
                    ax.legend()
                    plt.show()

            
            except Exception as e:
                print(f"<b>Error plotting {name}:</b> {e}")   
    
  
    hist2d_names = ["fHistptvsalpha", "fHistdEdx"] 
    for name in hist2d_names:
        if name in histograms:
            try:
               
                h2d = histograms[name]
                values, x_edges, y_edges = h2d.to_numpy()

             
                fig = plt.figure(figsize=(10, 8))
                gs = fig.add_gridspec(1, 1)
                ax = fig.add_subplot(gs[0, 0])

                
               
                hep.hist2dplot(values, x_edges, y_edges, ax=ax)
                hep.alice.label("Preliminary", data=True, lumi=50, com=13,rlabel="ALICE")    

                
                ax.set_title(name, fontsize=16)
                ax.set_xlabel("X Variable", fontsize=14)
                ax.set_ylabel("Y Variable", fontsize=14)
                plt.tight_layout()
                plt.savefig(f"{name}_2d.png", dpi=300, bbox_inches='tight')
                plt.show()
            except Exception as e:
                print(f"<b>Error plotting {name}:</b> {e}")


    plt.colorbar()
    plt.show()   

def main():
    
    histograms = process_root_file()
    
    if histograms:
    
        plot_histograms(histograms)

if __name__ == "__main__":
    main()

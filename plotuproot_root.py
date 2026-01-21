import sys
import uproot
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import m_e, c


def process_root_file(filename="coswindows.root"):
    try:
        # Open the ROOT file
        with uproot.open(filename) as file:
            # Print file contents
            print("<b>Keys in file:</b>")
            print(file.classnames())

            # Access the specific container
            container = file["MyTask/MyOutputContainer;1"]
            
            # List to store histograms of interest
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

                #other histos
            ]

            # Dictionary to store found histograms
            found_histograms = {}

            # Iterate through container to find histograms
            for obj in container:
                if hasattr(obj, "name") and obj.name in histograms_of_interest:
                    found_histograms[obj.name] = obj
                    print(f"<b>Found histogram:</b> {obj.name}")

            # Return the found histograms
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

    # Set ALICE style
    plt.style.use(hep.style.ALICE)

    # Mapping of histogram names to more descriptive titles
    histoname = {
        "cos0": "cos[-1;-0.5]",
        "cos1": "cos[-0.5,0]",
        "cos10": "cos[0;0.5]",
        "cos50": "cos[0.5,1]",
        "dEdx" : "Mult"
    }

    # Color palette for different histograms
    colors = ['blue', 'cyan', 'purple', 'magenta']

    # Create a single figure and GridSpec
    fig = plt.figure(figsize=(12, 8))  # Assign figure to variable
    gs = fig.add_gridspec(1, 1)        # Add GridSpec
    ax1 = fig.add_subplot(gs[0, 0])

  
    # Cos theta histis
    #canvas
    
    for (name, hist), color in zip([(k, v) for k, v in histograms.items() if k.startswith('cos')], colors):
        try:
        # Convert to histogram
            hist_converted = hist.to_hist()
            # Plot the histogram
            hep.histplot(hist_converted, 
                     ax=ax1, 
                     label=histoname.get(name, name), 
                     color=color)
        except Exception as e:
            print(f"<b>Error plotting {name}:</b> {e}")
       
     # Customize the plot
    ax1.set_title("InvMassLambdaP", fontsize=16)  # Use ax1.set_title()
    ax1.set_xlabel("Gev/c2", fontsize=14)        # Use ax1.set_xlabel()
    ax1.set_ylabel("Entries", fontsize=14)        # Use ax1.set_ylabel()
    # Add legend
    # Add grid
    ax1.grid(True, linestyle='--', alpha=0.5)
    # Adjust layout
    plt.tight_layout() 
    # Save the figure
    ax1.legend(title="CosTheta", loc='best')
    plt.savefig("cosine_windows_combined.png", dpi=300, bbox_inches='tight')   
    # Show the plot
    plt.show()
    #dedx 
  
    for name, hist in histograms.items():
        #other histos
        if not name.startswith('cos'):
            try:
                # Convert to histogram
                h = hist.to_hist()
                
                #  2d 
                if h.ndim == 2:
                    # Plot 2D histogram
                    fig = plt.figure(figsize=(10, 8))
                    ax = fig.add_subplot()
                    hep.hist2dplot(h, ax=ax)
                    
                    if name == "fHistdEdx":
                        # Define Bethe-Bloch function
                        def bethe_bloch(p, m):
                            beta = p / np.sqrt(p**2 + m**2)
                            gamma = 1 / np.sqrt(1 - beta**2)
                            return (1/beta**2) * (np.log(2*0.511*beta**2*gamma**2 / (23.2e-6)) - beta**2)
                        
                        #p = np.linspace(0.1, 10, 100)
                        ax.set_xscale('log')
                        # Ensure Bethe-Bloch curves cover the log range
                        p = np.logspace(np.log10(0.1), np.log10(10), 100)

                        for m, label, color in [(0.139, 'Pion', 'red'), (0.494, 'Kaon', 'blue'), (0.938, 'Proton', 'green')]:
                            ax.plot(p, bethe_bloch(p, m), color=color, label=label, linewidth=2)
                        
                    ax.legend()
                    ax.set_title(name)
                    plt.show()
                else:
                    # Plot 1D histogram
                    fig = plt.figure(figsize=(10, 7))
                    ax = fig.add_subplot()
                    hep.histplot(h, ax=ax, color='black')
                    ax.set_title(name)
                    ax.legend()
                    plt.show()
                    
            except Exception as e:
                print(f"<b>Error plotting {name}:</b> {e}")   
    
    # orIdentify 2D histograms (example names)
    hist2d_names = ["fHistptvsalpha", "fHistdEdx"]
    for name in hist2d_names:
        if name in histograms:
            try:
                # Load and convert 2D histogram
                h2d = histograms[name]
                values, x_edges, y_edges = h2d.to_numpy()

                # Create figure
                fig = plt.figure(figsize=(10, 8))
                gs = fig.add_gridspec(1, 1)
                ax = fig.add_subplot(gs[0, 0])

                # Plot using mplhep
                hep.hist2dplot(values, x_edges, y_edges, ax=ax)
                hep.alice.label("Preliminary", data=True, lumi=50, com=13,rlabel="ALICE")  # Use alice.label if available   
                # Customize
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
    # Process the ROOT file
    histograms = process_root_file()
    
    if histograms:
        # Plot the histograms
        plot_histograms(histograms)

if __name__ == "__main__":
    main()

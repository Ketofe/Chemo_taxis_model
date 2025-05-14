from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np



def animate(pos,nutrient_concentration, L, number_of_frames, interval=1):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Set axis limits
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)

    # Initial heatmap of the nutrient concentration
    im = ax.imshow(nutrient_concentration[0].T, extent=[0, L, 0, L],origin="lower")

    # Initial quiver plot for particle movement
    qv = ax.scatter(pos[0][:, 0], pos[0][:, 1], color="black",s=5)

    # Animation function
    def update(i):
        # Update nutrient concentration heatmap
        im.set_array(nutrient_concentration[i].T)

        # Update the position and direction of particles
        qv.set_offsets(pos[i])
 

        return im,          qv  
       # return qv
    anim = FuncAnimation(fig, update, frames=np.arange(1, number_of_frames), interval=interval, blit=False)

    plt.colorbar(im, label="Concentration")
  #  plt.show()
    
    return anim

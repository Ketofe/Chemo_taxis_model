import numpy as np
#from scipy.spatial import cKDTree
from findiff import FinDiff


def normalise_vector(vector):
        size=np.linalg.norm(vector)
        if size==0:
                  return vector
        else:
                  return vector/size
        




                                                         #Here one can set initial conditions if not they wil be set to default values further down
def Vicsek_with_communication_nutrient_diffusion3_5(v_0, N,L,c_D,evaporation,rate_of_addition, n_D,growth_rate,K,comsumption_rate,alpha,  Grid_size , iterations,initial_pos=None,initial_c_concentration=None,initial_n_concentration=None,print_iteration=True):
        orientation=np.zeros([iterations,N,2])
        pos=np.zeros([iterations,N,2])

        c_concentatrion=np.zeros([iterations,Grid_size,Grid_size])
        n_concentatrion=np.zeros([iterations,Grid_size,Grid_size])

        #For the gird
        xs=np.linspace(0,L,Grid_size)
        ys=np.linspace(0,L,Grid_size)
        
        spacing=xs[1]-xs[0]
        d2_dx2 = FinDiff(0, spacing, 2,periodic=True)  # d^2/dx^2
        d2_dy2 = FinDiff(1, spacing, 2,periodic=True)  # d^2/dy^2


        if initial_pos is None:
               initial_pos= [ [np.random.uniform(0,L),np.random.uniform(0,L) ] for k in range(N) ]
        pos[0]=initial_pos
       
           
                        

        if initial_c_concentration is not None:
                     c_concentatrion[0]=initial_c_concentration
              
        if initial_n_concentration is not None:
                     n_concentatrion[0]=initial_n_concentration

        

        for iteration in range(iterations-1):
                

                #Updating the grid with out the influence of the particles                 #The laplacian part                                                                 #The growth part                                           #the logistic part  
                c_concentatrion[iteration+1]=c_concentatrion[iteration]+c_D*(d2_dx2(c_concentatrion[iteration])+d2_dy2(c_concentatrion[iteration]))-evaporation*c_concentatrion[iteration] 
                n_concentatrion[iteration+1]=n_concentatrion[iteration]+n_D*(d2_dx2(n_concentatrion[iteration])+d2_dy2(n_concentatrion[iteration]))+growth_rate*n_concentatrion[iteration]*(1-K*n_concentatrion[iteration])

             
               
            
   
                
                #Now looping over all the particles and doing the angle update
                for particle_index in range(N):
                      #Now findig the grid idexes of the particle
                      x_index=int(pos[iteration][particle_index][0]//spacing)
                      y_index=int(pos[iteration][particle_index][1]//spacing)

                      #Now updating the concenctracions

                      c_concentatrion[iteration+1][x_index][y_index]+=rate_of_addition 
                                                                      #Using plus here since it is also updated outside particle loop
                      n_concentatrion[iteration+1][x_index][y_index]= n_concentatrion[iteration+1][x_index][y_index]-comsumption_rate*c_concentatrion[iteration][x_index][y_index] 


                      #Now updating the oreintation based on the gradient

                   #   c_gradient_x=c_concentatrion[iteration][(x_index+1)%Grid_size][y_index]-c_concentatrion[iteration][x_index][y_index] 
                   #   c_gradient_y=c_concentatrion[iteration][x_index][(y_index+1)%Grid_size]-c_concentatrion[iteration][x_index][y_index]
                   #   n_gradient_x=n_concentatrion[iteration][(x_index+1)%Grid_size][y_index]-n_concentatrion[iteration][x_index][y_index]
                   #   n_gradient_y=n_concentatrion[iteration][x_index][(y_index+1)%Grid_size]-n_concentatrion[iteration][x_index][y_index]
                      
                      #Central derivative
                      c_gradient_x = (c_concentatrion[iteration][(x_index+1)%Grid_size][y_index] - 
                           c_concentatrion[iteration][(x_index-1)%Grid_size][y_index]) 
                      c_gradient_y = (c_concentatrion[iteration][x_index][(y_index+1)%Grid_size] - 
                                   c_concentatrion[iteration][x_index][(y_index-1)%Grid_size])
                      n_gradient_x = (n_concentatrion[iteration][(x_index+1)%Grid_size][y_index] - 
                           n_concentatrion[iteration][(x_index-1)%Grid_size][y_index]) 
                      n_gradient_y = (n_concentatrion[iteration][x_index][(y_index+1)%Grid_size] - 
                                   n_concentatrion[iteration][x_index][(y_index-1)%Grid_size]) 

                       
                      c_gradient=np.array([c_gradient_x,c_gradient_y])
                      n_gradient=np.array([n_gradient_x,n_gradient_y])
                                                              #Normalising to ensuree that it has a size of 1 or 0. The sum will not have a size of 1 or 0 if one of the gradients is 0
                      orientation[iteration][particle_index]=normalise_vector((1-alpha)*normalise_vector(c_gradient)+alpha*normalise_vector(n_gradient))

                if print_iteration: 
                        print('iteration:',iteration)

                #Also doing this if statement for orientation
               
                #Updating the postions 
                pos[iteration + 1] = (pos[iteration] + v_0 * orientation[iteration]) % L
        
        #Now updatin only orientations since orientation get update like iteration and not iteration+1 so using the last iteration here
        iteration=-1
        for particle_index in range(N):
                      #Now findig the grid idexes of the particle
                      x_index=int(pos[iteration][particle_index][0]//spacing)
                      y_index=int(pos[iteration][particle_index][1]//spacing)

                       #Central derivative
                      c_gradient_x = (c_concentatrion[iteration][(x_index+1)%Grid_size][y_index] - 
                           c_concentatrion[iteration][(x_index-1)%Grid_size][y_index]) 
                      c_gradient_y = (c_concentatrion[iteration][x_index][(y_index+1)%Grid_size] - 
                                   c_concentatrion[iteration][x_index][(y_index-1)%Grid_size])
                      n_gradient_x = (n_concentatrion[iteration][(x_index+1)%Grid_size][y_index] - 
                           n_concentatrion[iteration][(x_index-1)%Grid_size][y_index]) 
                      n_gradient_y = (n_concentatrion[iteration][x_index][(y_index+1)%Grid_size] - 
                                   n_concentatrion[iteration][x_index][(y_index-1)%Grid_size]) 

                       
                      c_gradient=np.array([c_gradient_x,c_gradient_y])
                      n_gradient=np.array([n_gradient_x,n_gradient_y])
                                                              #Normalising to ensuree that it has a size of 1 or 0. The sum will not have a size of 1 or 0 if one of the gradients is 0
                      orientation[iteration][particle_index]=normalise_vector((1-alpha)*normalise_vector(c_gradient)+alpha*normalise_vector(n_gradient))


                      
        return pos,orientation,c_concentatrion,n_concentatrion



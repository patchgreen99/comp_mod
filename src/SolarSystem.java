/**
 * SolarSystem class to encapsulate the physical interactions between objects Particle3D
 * Class Variables :  private Particle3D[] particle_array, private double total_energy, private static double total_initial_momentum, private static double com_Velocity
 */
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;

public class SolarSystem {

/**
 * Object SolarSystem Elements
 */
private Particle3D[] particleArray;
private double totalEnergy; 
private static double totalInitialMomentum;
private static double comVelocity;
private int numberTimesteps;
private double timestep;

/**
 * Constructor for SolarSystem which contructs an array of particles read from a file and adjusts intial velocities
 * @param  file a String which is the name of the file containing the data on the particles
 * @return  the object
 */
	public SolarSystem(String file){

	}

/**
 * Sets the particle at index particleIndex to the argument particle
 * @param particleIndex is an integer 
 * @param particle      the particle3D to set it as
 */
	public void setParticle(int particleIndex, Particle3D particle){
		particleArray[particleIndex] = particle;
	}

/**
 * Sets the total energy of the particle
 * @param energy gives a double and is the energy of the particle
 */
	public void setTotalEnergy(double energy){
		totalEnergy = energy;
	}

/**
 * Gets the Particle3D at index particleIndex
 * @param  particleIndex  the index of the Particle3D to return
 * @return Particle3D     the particle3D in the particle array with index particleIndex
 */
	public Particle3D getParticle(int particleIndex){
		return particleArray[particleIndex];	
	}

/**
 * Gets the total energy of the SolarSystem
 * @return a double of the total energy
 */
	public double getTotalEnergy(){
		return totalEnergy;
	}

/**
 * Subtract the centre of mass velocity from each particles velocity so that the simulation is in the frame of reference of space and not the sun (only at the start)
 */
	public void com_Correction(){

	}

/**
 * Creates the string to Write the position of each particle and details of which particle it is for each timestep. It returns 2 /n Point = 1 /n s11 x11 y11 z11 where 2 is the index of the plot point is the index of timesteps and the tuple geolocates the particle
 * @return 2 /n Point = 1 /n s11 x11 y11 z11 where 2 is the index of the plot point is the index of timesteps and the tuple geolocates the particle
 */
	public String toString(){
		return "";
	}


/**
 * Returns the seperation vector of the orbit for particle with index particle 
 * @param  integer equal to the particle index
 * @return returns the orbital radius of that particle
 */
	public Vector3D getOrbitalSeperation(int particle){
		return new Vector3D();
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * updates the velocity of the particle with new timstep
 * @param timestep is a double of the updated timestep
 */
	public void updateVelocity(double timestep){
		
	}

/**
 * updates the position of the particle given the new timestep
 * @param timestep is a double of which to update the timestep
 */
	public  void updatePosition(double timestep){

	} 

/**
 * Updates the force of each particle given the new velocity, position and sums the gravitational attractive force from every other particle in the SolarSystem
 */
	public void updateForce(){

	}


/**
 * updates the energy of the SolarSystem given the new particle positions and velocities
 */
	public void updateEnergy(){

	}

/**
 * Returns the half time period taken by the particle with index particle to orbit the star in the SolarSystem. This is done by recording the time at every minimum arccosine of the normalised dot product of the position of the star and particle. 
 * @param  integer of the particle in question
 * @param  fractionOfOrbit fraction of the orbit
 * @return returns the number of orbits
 */
	public int timePeriod(int particle){
		return 1;
	}

/**
 * Gives the aphelion (point where the particle is farthest from the Sun) saves the seperation vector magnitude for particle with index particle per itteration if larger than the value previously saved
 * @param  particle integer of the particle
 * @return the aphelion
 */
	public double aphelion(int particle){
		return 1.0;
	}

/**
 * Gives the perhelion (point where the particle is closest to the Sun) saves the seperation vector magnitude for particle with index particle per itteration if smaller than the value previously saved
 * @param  particle integer of the particle
 * @return the perihelion
 */
	public double perihelion(int particle){
		return 1.0;
	}



	/**
	 * Constructor that sets up an array of Particle3D objects taken from ( SUN input file what ever we want ,  can have v = 0)
	 * SolarSystem class to store data for particles in the system and then simulates the system over time using integration algorithms and writes a trajectory file with the simulation results so that these can be visualised using VMD
	 * @param  a         object SolarSystem in which the simulation will take place in
	 * @param  num_steps integer that gives number of steps in integration
	 * @param  step      double that gives the magnitude of each timestep in the integration
	 * @return           returns object N_body_simulation
	 */
		public N_body_simulation(SolarSystem a , int num_steps, double step){

		}

	/**
	 * Contructs a simulation that is the gravitational interaction of N bodies by reading the planetary data and simulation constants file.
	 * @param  sim_file    name of a file in which the simulation parameters are written to once they have been calculated
	 * @param  planet_File name of a file containing the details of the bodies in the system
	 * @return             a simulation of the solar system as a set of parameters that can be visualised using VMD
	 */
		public N_body_simulation(String sim_file, String planet_File){

		}

	/**
	 * returns integer number of timesteps set in the simulation is retrieved
	 * @return integer representing the number of timesteps
	 */
		public int getNumberTimesteps(){
			return 1;
		}

	/**
	 * returns the double timestep between each iteration of the integration is retrieved from the simulation
	 * @return double representing the time step
	 */
		public double getTimestep(){
			return 1.0;
		}

	/**
	 * A trajectory file is written that contains the trajectory calculated using the input parameters
	 * @param filename name of the file created. read from terminal.
	 */
		public void outputTrajectoryFile(String filename){

		}

	/**
	 * file is created to store the calculated fluctuations in energy for each body in the system
	 * @param filename name of the file created. read from terminal.
	 */
		public void outputEnergyFluctuations(String filename){

		}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * A verlet update in which the velocity and position of each planet is updated per time step
	 * @param timestep double equal to the timestep of the simulation
	 */

	    public static void verletUpdate(double timestep){
		
		// Update the position using current velocity and force
		//orbiter.updatePosition(timestep,TimeIntegral.getForce());
		// Update forces based on new positions
		//Vector3D force_new = TimeIntegral.getGravity(orbiter,orbitee);
		
		// Update the velocity, based on average of current and new force
		//orbiter.updateVelocity(timestep,Vector3D.addVector3D(TimeIntegral.getForce(), force_new).scalarMul(0.5));
		
		// Update the force variable to the new value
		//TimeIntegral.setForce(force_new);
	    }
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

	/**
	 * Main method. Takes 3 filenames from the terminal as an argument. Planetary Data , Simulation Constants , Output file
	 * @param  args        file names with particle data
	 * @throws IOException ensures that program goes into a try/catch block where program can be checked for errors and rectified
	 */
		public static void main(String[] args) throws IOException {

		}
	}

}





























 

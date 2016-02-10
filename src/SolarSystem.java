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
private Particle3D[] particle_array;
private double total_energy; 
private static double total_initial_momentum;
private static double com_Velocity;

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

	}

/**
 * Sets the total energy of the particle
 * @param energy gives a double and is the energy of the particle
 */
	public void setTotalEnergy(double energy){

	}

/**
 * Gets the Particle3D at index particleIndex
 * @param  particleIndex  the index of the Particle3D to return
 * @return Particle3D     the particle3D in the particle array with index particleIndex
 */
	public Particle3D getParticle(int particleIndex){
		return new Particle3D();
	}

/**
 * Gets the total energy of the SolarSystem
 * @return a double of the total energy
 */
	public double getTotalEnergy(){
		return 1.0;
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



}





























 

/**
 * N_body_simulation encapsulates the simulation in time a object SolarSystem
 * Class Variables:  private SolarSystem N_bodies, private int numberTimesteps, private double timestep;
 */

import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;

public class N_body_simulation {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Object N_body_simulation Elements
 */

private SolarSystem N_bodies;
private int numberTimesteps;
private double timestep;

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

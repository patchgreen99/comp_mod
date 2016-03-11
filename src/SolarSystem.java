/**
 * SolarSystem class to encapsulate the physical interactions between objects Particle3D
 * Class Variables :  private Particle3D[] particle_array, private double total_energy, private static double total_initial_momentum, private static double com_Velocity
 */
import java.util.ArrayList;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.Math;

public class SolarSystem {

	/**
	 * Object SolarSystem Elements
	 */
	private Particle3D[] particleArray;
	private double totalEnergy; 
	private int planets;
	private double bigG;
	private ArrayList<String> labels;
	private double timestep;
	

	/**
	 * Constructor for SolarSystem which constructs an array of particles read from a file and adjusts initial velocities
	 * @param  file a String which is the name of the file containing the data on the particles
	 * @return  the object
	 * @throws FileNotFoundException 
	 */
	public SolarSystem(String fileName, double t) throws FileNotFoundException{
		timestep = t;
		labels = new ArrayList<String>();
		BufferedReader file = new BufferedReader(new FileReader(fileName));
		Scanner scan = new Scanner(file);
		planets = scan.nextInt();
		particleArray = new Particle3D[planets];
		
		for( int i = 0 ; i < planets ; i++){
			this.particleArray[i] = new Particle3D(scan);
		}
		for( int i = 0 ; i < planets ; i++){
			labels.add(this.particleArray[i].getLabel());
		}
		scan.close();
		correctForCOM();
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
	 *
	 * @param 
	 */
	public void setG(double g){
		bigG = g;
	}

	/**
	 * 
	 * @return
	 */
	public double getG(){
		return bigG;
	}

	/**
	 * Gets the Particle3D at index particleIndex
	 * @param  particleIndex  the index of the Particle3D to return
	 * @return Particle3D     the particle3D in the particle array with index particleIndex
	 */
	public Particle3D getParticle(String label){
		int i = labels.indexOf(label);
		return particleArray[i];	
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
	public void correctForCOM(){
		Vector3D momentumSum = new Vector3D();
		double massSum = 0;
		for(int i = 0 ; i < particleArray.length ; i++){
			massSum = massSum + particleArray[i].getMass();
			momentumSum = Vector3D.addVector3D(momentumSum, particleArray[i].getVelocity().scalarMul(particleArray[i].getMass()));
		}
		
		Vector3D velocityCorrection = momentumSum.scalarDiv(massSum);
		
		for(int j = 0 ; j < particleArray.length ; j++){
			Vector3D temp = particleArray[j].getVelocity();
			particleArray[j].setVelocity(Vector3D.subVector3D(temp, velocityCorrection));
		}
	}

	/**
	 * Creates the string to Write the position of each particle and details of which particle it is for each timestep. It returns 2 /n Point = 1 /n s11 x11 y11 z11 where 2 is the index of the plot point is the index of timesteps and the tuple geolocates the particle
	 * @return 2 /n Point = 1 /n s11 x11 y11 z11 where 2 is the index of the plot point is the index of timesteps and the tuple geolocates the particle
	 */
	public String toString(){
		String answer = "";
		for(int i = 0 ; i < particleArray.length ; i++){
			answer = answer + particleArray[i].toString()+"\n" ;
		}
		return answer ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * Updates the force of each particle given the new velocity, position and sums the gravitational attractive force from every other particle in the SolarSystem
	 */
	public void updateAllForce(){
		Vector3D[] forceSumArray = new Vector3D[particleArray.length];
		for(int k = 0 ; k < particleArray.length ; k++){
			forceSumArray[k] = new Vector3D();
		}
		
		for(int i = 0 ; i < particleArray.length-1 ; i++){	
			for(int j = i+1 ; j < particleArray.length ; j++){
				Vector3D tempForce = Particle3D.getGravitationalAttraction(particleArray[i],particleArray[j],this.getG());
				Vector3D iForce = Vector3D.addVector3D(forceSumArray[i], tempForce);
				Vector3D jForce = Vector3D.addVector3D(forceSumArray[j],tempForce.scalarMul(-1.0));
				forceSumArray[i] = iForce;
				forceSumArray[j] = jForce;
			}
		}
		
		for(int l = 0 ; l < particleArray.length ; l++){
			particleArray[l].setForce(forceSumArray[l]);
		}
	}
	
	/**
	 * 
	 */
	public void updatePosition(double ts){
		for(int i = 0 ; i < particleArray.length ; i++){
			particleArray[i].updatePosition(ts);
		}
	}
	
	/**
	 * Updates the force of each particle given the new velocity, position and sums the gravitational attractive force from every other particle in the SolarSystem
	 */
	public void updateVelocity(double ts){
		// Update the position using current velocity and force
		// Update forces based on new positions
		// Update the velocity, based on average of current and new force
		// Update the force variable to the new value
		
		for(int i = 0 ; i < particleArray.length; i++){	
			//particleArray[i].updatePosition(ts);
			Vector3D forceSum = new Vector3D();
			
			for(int j = 0 ; j < particleArray.length ; j++){
				if(j != i){
					Vector3D newForce = Particle3D.getGravitationalAttraction(particleArray[i],particleArray[j],this.bigG);
					Vector3D tempForce = Vector3D.addVector3D(forceSum, newForce);
					forceSum = tempForce;
				}
			}			
			particleArray[i].updateVelocity(ts,Vector3D.addVector3D(particleArray[i].getForce(), forceSum).scalarMul(0.5));
			particleArray[i].setForce(forceSum);
		}

	}


	/**
	 * updates the energy of the SolarSystem given the new particle positions and velocities
	 */
	public void updateEnergy(){
		double sumEnergy = 0;
		for(int k = 0 ; k < particleArray.length ; k++){
			sumEnergy = sumEnergy + particleArray[k].getKE();
		}
		for(int i = 0 ; i < particleArray.length-1 ; i++){
			for(int j = i + 1 ; j < particleArray.length ; j++){
				sumEnergy = sumEnergy + Particle3D.getGravEnergy(particleArray[i],particleArray[j],this.getG());
			}
		}
		totalEnergy = sumEnergy;
	}


	/**
	 * A trajectory file is written that contains the trajectory calculated using the input parameters
	 * @param filename name of the file created. read from terminal.
	 * @throws IOException 
	 */
	public void outputTrajectoryFile(PrintWriter output, int i) throws IOException{
	    output.print(planets+"\n"); // colour code
	    output.print("Point = "+i+"\n");
	    output.print(toString());
	}

	/**
	 * file is created to store the calculated fluctuations in energy for each body in the system
	 * @param filename name of the file created. read from terminal.
	 */
	public void outputEnergyFluctuations(PrintWriter output, int i){
		output.print(i+"	"+this.getTotalEnergy()+"\n");
	}


	public static double getOrbitalAngle(Particle3D a, Particle3D b){
		if (a.getMass() > b.getMass()){
			return Particle3D.orbitAngle(b, a);
		}else{
			return Particle3D.orbitAngle(a, b);
		}
	}
	
	public void updateOrbitData(String a, String b, int i){
		double currentAngle = SolarSystem.getOrbitalAngle(getParticle(a), getParticle(b));
		Particle3D orbiting;
		Particle3D star;
		if (getParticle(a).getMass() < getParticle(b).getMass()){
			orbiting = getParticle(a);
			star = getParticle(b);
		}else{
			orbiting = getParticle(b);
			star = getParticle(a);
		}
		
		//ORBIT CALCULATION
		//////////////////////////////////////////////////////
		orbiting.addAngleStack(currentAngle);
		if (orbiting.getAngleStack().size() > 3){
			orbiting.removeAngleStack();
			if (orbiting.getAngleStack().get(1) < orbiting.getAngleStack().get(0) && orbiting.getAngleStack().get(1) < orbiting.getAngleStack().get(2) ){ 
				orbiting.addTimeStack((i-1)*timestep);
				if (orbiting.getTimeStack().size() > 2){
					orbiting.removeTimeStack();
					if (orbiting.getTimeStack().get(1) - orbiting.getTimeStack().get(1) < orbiting.getPeriod()){
						orbiting.setPeriod(orbiting.getTimeStack().get(1) - orbiting.getTimeStack().get(0));
					}
				}
			}
		}
		
		//DISTANCE CALCULATIONS
		//////////////////////////////////////////////
		if(orbiting.getAphelion() < Particle3D.seperationAway(orbiting,star).magnitude()){
			orbiting.setAphelion(Particle3D.seperationAway(orbiting,star).magnitude());
		}
		
		if(orbiting.getPerihelion() > Particle3D.seperationAway(orbiting,star).magnitude()){
			orbiting.setPerihelion(Particle3D.seperationAway(orbiting,star).magnitude());
		}
	}
	

	/**
	 * Main method. Takes 3 filenames from the terminal as an argument. Planetary Data , Simulation Constants , Output file
	 * @param  args        file names with particle data
	 * @throws IOException ensures that program goes into a try/catch block where program can be checked for errors and rectified
	 */
	public static void main(String[] args) throws IOException {
		BufferedReader file = new BufferedReader(new FileReader("param.dat"));
		Scanner scan = new Scanner(file);
		double timestep = scan.nextDouble();
		int numTimeStep = scan.nextInt();
		double bigG = scan.nextDouble();
		scan.close();
		
		SolarSystem s = new SolarSystem("input.dat", timestep);
		s.setG(bigG);
		
		PrintWriter outputTrajectory = new PrintWriter(new FileWriter("outputTrajectory.xyz"));
		PrintWriter outputEnergy = new PrintWriter(new FileWriter("outputEnergy.dat"));
		
		////////////////////////////////////////////////////
		//Movement
		////////////////////////////////////////////////////
		s.updateAllForce();
		
		for(int i=0 ; i < numTimeStep ; i++){
			s.outputTrajectoryFile(outputTrajectory, i);
			s.outputEnergyFluctuations(outputEnergy, i);
			
			s.updatePosition(timestep);
			s.updateVelocity(timestep);
			s.updateEnergy();
			
			s.updateOrbitData("SUN", "EARTH", i);
			s.updateOrbitData("EARTH", "MOON", i);
					
		}
		
		outputTrajectory.close();
		outputEnergy.close();
		System.out.print("Earth has period "+s.getParticle("EARTH").getPeriod()/(86400)+" in days" + "\n");
		System.out.print("Earth has Aphelion "+s.getParticle("EARTH").getAphelion()+" in meters" + "\n");
		System.out.print("Earth has Perihelion "+s.getParticle("EARTH").getPerihelion()+" in meters" + "\n");
		System.out.print("Earth has Eccentricity "+s.getParticle("EARTH").getEccentricity()+ "\n");
		System.out.print("Earth has Semi Major Axis "+s.getParticle("EARTH").getSemiMajorAxis()+ " in meters" + "\n");
		
		//System.out.print("Moon has period "+s.getParticle("MOON").getPeriod()/(86400)+" in days" + "\n");
	}


}






























 

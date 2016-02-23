import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;


public class TimeIntegralVerlet{

//Class variables bigG, force and total Energy which are specific to the simulation 
	
    private static double bigG = 1.0;
    private Vector3D force;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Update the force for that simulation
    public void setForce(Vector3D f){
    	this.force = f;}

    //obtain the force for that simulation
    public Vector3D getForce(){
    	return this.force;}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Constructor for the time integral.
    public TimeIntegralVerlet(){}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Simulation outputs spacial coordinates and Energy against time (x,y) (t,E)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static void main(String[] args) throws IOException {
		
    // Open the input file
    String fileName = "input.dat";
    BufferedReader file = new BufferedReader(new FileReader(fileName));
    Scanner scan = new Scanner(file);
    	
	// Open the output file for Position
    String outFilePosition = "outputPositionVerlet.dat";
    PrintWriter outputPosition = new PrintWriter(new FileWriter(outFilePosition));
    
    // Open the output file for Energy
    String outFileEnergy = "outputEnergyVerlet.dat";
    PrintWriter outputEnergy = new PrintWriter(new FileWriter(outFileEnergy));
   
    //Construct a simulation
    TimeIntegralVerlet TimeIntegral = new TimeIntegralVerlet();

    //construct two different particles	from a file
	Particle3D orbiter = new Particle3D(scan);
	Particle3D orbitee = new Particle3D(scan);
	scan.close();
		
	// set the initial force between the two particles
	TimeIntegral.setForce(TimeIntegral.getGravity(orbiter,orbitee));

////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
	
		// Number of time steps
		double numstep = 100;
		// Size of time step
		double timestep = 0.1;
		// Initial time
		double t = 0;	

		// Print the initial data about the particle to file
		outputPosition.printf("%10.5f %10.5f\n",orbiter.getx().getX(),orbiter.getx().getY());
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Start Loop over time
		for (int i=0;i<numstep;i++){
	    	    	
			TimeIntegral.verletUpdate(TimeIntegral,orbiter,orbitee,timestep);
	    	
			// Increase the time
			t = t + timestep;
	    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Print the current information on particle to file
			outputPosition.printf("%10.5f %10.5f\n",orbiter.getx().getX(),orbiter.getx().getY());
			outputEnergy.printf("%10.5f %10.8e\n",t,TimeIntegral.totEnergy(orbiter, orbitee));
		//end loop over time
		}
	// Close the output file
	outputPosition.close();
	outputEnergy.close();
    }//end of the simulation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //TimeIntegral METHODS
    // Obtain gravitational force
    // Obtain total energy
    // run a Verlet update
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static Vector3D getGravity(Particle3D orbForced, Particle3D orbForcing){
	return Vector3D.subVector3D(orbForced.getx(), orbForcing.getx()).scalarMul
	    (bigG*orbForced.getMass()*orbForcing.getMass()*-1.0*
	     Math.pow(Vector3D.subVector3D(orbForced.getx(), orbForcing.getx()).magnitude(),-3.0));
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static double totEnergy(Particle3D orbit, Particle3D station){
    	return orbit.getKE() - bigG*orbit.getMass()*station.getMass()
	    /Vector3D.subVector3D(orbit.getx(), station.getx()).magnitude();
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static void verletUpdate(TimeIntegralVerlet TimeIntegral,Particle3D orbiter, Particle3D orbitee,double timestep){
	
	// Update the position using current velocity and force
	orbiter.updatePosition(timestep,TimeIntegral.getForce());
	// Update forces based on new positions
	Vector3D force_new = TimeIntegral.getGravity(orbiter,orbitee);
	
	// Update the velocity, based on average of current and new force
	orbiter.updateVelocity(timestep,Vector3D.addVector3D(TimeIntegral.getForce(), force_new).scalarMul(0.5));
	
	// Update the force variable to the new value
	TimeIntegral.setForce(force_new);
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

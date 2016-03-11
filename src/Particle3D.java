/**
 * Particle3D class to obtain the energy of the particle, its force and distance between particles 
 * Class variables :     private Vector3D v, private Vector3D x, private double mass, private String label, private Vector3D force
 */
import java.util.ArrayList;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;
import java.math.BigDecimal;
import java.text.DecimalFormat;

//Using Vector3D to define a particle's velocity, mass and position in vector components.
public class Particle3D {
    
    /*
     * Object vector elements
     */
  
    private Vector3D v;
    private Vector3D x;
    private double mass;
    private String label;
    private Vector3D force;
    private ArrayList<Double> angleStack;
	private ArrayList<Double> timeStack;
	private double period;
	private double aphelion;
	private double perihelion;
    
    /* 
     * Default constructor to set the Particle3D to zero
     */
//Default 3D particle
    public Particle3D(){
        x = new Vector3D();
        v = new Vector3D();
        mass = 0;
        label = null;
        force = new Vector3D();
        angleStack = new ArrayList<Double>();
		timeStack = new ArrayList<Double>();
		period = 1000000000;
		aphelion = 0;
		perihelion = 10E100;
    }
    
    /**Sets the velocity vector of the particle
     * 
     * @param v is a 3D velocity vector and sets the velocity vector
     */
    
    public void setVelocity(Vector3D v){
        this.v = v;
    }
    
    /**Sets the position vector of the particle.
     * 
     * @param x is a 3D position vector and set the position vector
     */
    public void setPosition(Vector3D x){
        this.x = x;
    }
    
    /**Sets the mass of the particle
     * 
     * @param m is a double and sets the mass of the particle
     */
    
    
    public void setMass(double m){
        this.mass = m;
    }
    
    /**Sets the label for the particle
     * 
     * @param s is a String and set the label of the particle
     */

    public void setLabel(String s){
        this.label = s;
    }

    /**
     * sets the force of object Particle3D
     * @param f the new force
     */
    
    public void setForce(Vector3D f){
        this.force = f;
    }
    
    public void addTimeStack(Double e){
        this.timeStack.add(e);
    }
    
    public void addAngleStack(Double e){
        this.angleStack.add(e);
    }
    
    public void removeTimeStack(){
        this.timeStack.remove(0);
    }
    
    public void removeAngleStack(){
        this.angleStack.remove(0);
    }
    
    public void setPeriod(double e){
        this.period = e;
    }
    
    public void setAphelion(double e){
        this.aphelion = e;
    }
    
    public void setPerihelion(double e){
        this.perihelion = e;
    }
    
    
    ///////////////////////////////////////////////////////
    //Getter method for velocity, position, mass and the label of the 3D particle
    
    /**Gets the velocity vector of Particle3D 
     * 
     * @return a velocity vector instance representing the Particle3D velocity Vector.
     */
    
    public Vector3D getVelocity(){
        return this.v;
    }
    
    /**Gets the position vector of Particle3D 
     * 
     * @return a position vector instance representing the Particle3D position Vector.
     */
    
    public Vector3D getPosition(){
        return this.x;
    }
    
    /**Gets the mass of Particle3D
     * 
     * @return double instance representing the mass of Particle3D.
     */
    
    public double getMass(){
        return this.mass;
    }
    
    /**Gets the label of Particle3D
     * 
     * @return String of the Particle3D
     */
    
    public String getLabel(){
        return this.label;
    }

    /**
     * gets the force for Partice3D
     * @return Vector3D the force
     */
    public Vector3D getForce(){
        return this.force;
    }
    
    public ArrayList<Double> getTimeStack(){
        return this.timeStack;
    }
    
    public ArrayList<Double> getAngleStack(){
        return this.angleStack;
    }
    
    public double getPeriod(){
        return this.period;
    }
    
    public double getAphelion(){
        return this.aphelion;
    }
    
    public double getPerihelion(){
        return this.perihelion;
    }
    
    public double getEccentricity(){
        return (this.aphelion-this.perihelion)/(this.aphelion+this.perihelion);
    }
    
    public double getSemiMajorAxis(){
        return getAphelion()/(1+getEccentricity());
    }
    /////////////////////////////////////////////////////////////
    
    //Prints particles position in a readable form
    
    /**Returns a String which represents the position of Particle3D.
     * 
     */
    public String toString(){
        String tempLabel = this.getLabel();
        double x1_str = this.getPosition().getX();
        double x2_str = this.getPosition().getY();
        double x3_str = this.getPosition().getZ();
        DecimalFormat format = new DecimalFormat("0.00E00");
        return tempLabel+" "+format.format(x1_str)+" "+format.format(x2_str)+" "+format.format(x3_str);
    }
    /////////////////////////////////////////////////////////////

    //Scans values (Label, Mass, x, y, z, xvelocity, yvelocity, zvelocity) for components from a file and prints them
    
    /**Scans the values for Particle 3D (the label, mass, position vector and velocity vector)
     * 
     * @param scan gives the values scanned from the file.
     */
    public Particle3D(Scanner scan){
        label = (String) scan.next();
        mass = scan.nextDouble();
        x = new Vector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble()).scalarMul(1.49597870700E+11); //meters
        v = new Vector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble()).scalarMul(1.49597870700E+11/86400.0); // meters per second
        force = new Vector3D();
        angleStack = new ArrayList<Double>();
		timeStack = new ArrayList<Double>();
		period = 1000000000;
		aphelion = 0;
		perihelion = 10E100;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Finding the particles kinetic energy
    
   /**Gets the kinetic energy of the particles
    * 
    * @return a double representing the kinetic energy of the particle using it's mass, 
    * position vector and velocity vector.
    */
    public double getKE(){
        return 0.5*this.getMass()*(this.getVelocity().getX()*this.getVelocity().getX()+this.getVelocity().getY()*this.getVelocity().getY()
                +this.getVelocity().getZ()*this.getVelocity().getZ());
    }
    //updating the velocity knowing the force and timestep
    
    /**Updates the velocity using known force and timestep.
     * 
     * @param timestep is a double which updates timestep
     * @param force is a vector which updates the force
     */
    public void updateVelocity(double timestep, Vector3D tempForce){
    	Vector3D tempVelocity = this.getVelocity();
        this.setVelocity(Vector3D.addVector3D(tempVelocity, tempForce.scalarMul(timestep/this.getMass())));
    }
    //update position knowing the force and timestep
    /**Updates the position given the force and the timestep
     * 
     * @param timestep is a double that updates the position
     * @param force is a vector that updates the position
     */
    public void updatePosition(double timestep){
    	Vector3D tempPosition = this.getPosition();
        this.setPosition(Vector3D.addVector3D(Vector3D.addVector3D(tempPosition,this.getVelocity().scalarMul(timestep)),this.getForce().scalarMul(0.5*timestep*timestep/this.getMass())));
    }
    //distance between two particles   
    /**Calculates the seperation vector between two particles VECTOR FROM b TO a
     * 
     * @param a is the reference Particle3D
     * @param b is the Particle3D a distance away
     * @return the Vector3D  between the two particles 
     */
    public static Vector3D seperationAway(Particle3D a, Particle3D b){
        return Vector3D.subVector3D(a.getPosition(),b.getPosition());
    }

    //Static Methods

    /**
     * gravitational force on the particle with respect to the star pointing from moon to star
     * @param  particle Particle3D on which the force is acting 
     * @param  star   Particle3D by which the attraction is happening
     * @param  bigG   [description]
     * @return        Vector3D the force that is acting on the particle
     */
    public static Vector3D getGravitationalAttraction(Particle3D moon, Particle3D star, double bigG){
        return Particle3D.seperationAway(moon, star).scalarMul
        	    (bigG*moon.getMass()*star.getMass()*-1.0*
        	   	     Math.pow(Particle3D.seperationAway(moon, star).magnitude(),-3.0));
    }


    /**
     * The total energy for the 2 particle system
     * @param  particle one Particle3D in the system
     * @param  star   another Particle3D in the system
     * @param  bigG   [description]
     * @return        a double equal to the total energy for that pair of particles
     */
    public static double getGravEnergy(Particle3D planet1, Particle3D planet2, double bigG){
    	return  - (bigG*planet1.getMass()*planet2.getMass())
    		    /(Particle3D.seperationAway(planet1, planet2).magnitude());   
    }
    
    /**
     * Returns the angle of orbit from -180 to 180 with zero on the sun's position vector 
     * @param  integer of the particle in question
     * @param  fractionOfOrbit fraction of the orbit
     * @return returns the number of orbits
     */
    	public static double orbitAngle(Particle3D moon, Particle3D sun){
    		Vector3D seperation = Particle3D.seperationAway(sun, moon); //sun to moon
    		double cosine = Vector3D.dotVector3D(seperation,sun.getPosition().scalarMul(-1.0))/(seperation.magnitude()*sun.getPosition().magnitude());
    		double sine = Vector3D.crossVector3D(seperation,sun.getPosition().scalarMul(-1.0)).magnitude()/(seperation.magnitude()*sun.getPosition().magnitude());
    		if (cosine >= 0){
    				//System.out.print(Math.toDegrees(Math.asin(sine))+"\n");
                    return Math.toDegrees(Math.asin(sine));
    		}else{
                    double temp_angle = Math.toDegrees(Math.asin(sine));
                    if (temp_angle >= 0){
                    	//System.out.print(180 - temp_angle+"\n");
                        return 180 - temp_angle;
                    }
                    else{
                    	//System.out.print(-(180 + temp_angle)+"\n");
                        return -(180 + temp_angle);
                    }
    		}
    	}


    /**
     * sets velocity of moon to do circular orbit
     * @param  particle the orbiting particle with the outputted velocity
     * @param  star   the particle its orbiting around 
     */
    public static void setCircularOrbit(Particle3D moon, Particle3D star){
    	Vector3D moonDirection = moon.getVelocity().scalarDiv(moon.getVelocity().magnitude());
    	double factor = Math.sqrt((moon.getMass()+star.getMass())/Particle3D.seperationAway(moon, star).magnitude());
    	moon.setVelocity(moonDirection.scalarMul(factor)); 
    }    


}



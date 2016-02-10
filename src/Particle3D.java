/**
 * Particle3D class to obtain the energy of the particle, its force and distance between particles 
 * Class variables :     private Vector3D v, private Vector3D x, private double mass, private String label, private Vector3D force
 */
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;

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
    
    /* 
     * Default constructor to set the Particle3D to zero
     */
//Default 3D particle
    public Particle3D(){
        this.x = new Vector3D();
        this.v = new Vector3D();
        this.mass = 0;
        this.label = null;
    }
    /**Explicit constructor. Constructs a new Particle3D from given x, y and z components, it's mass and label
     * 
     * @param X is a 3D vector of it's position.
     * @param V is a 3D vector of its velocity.
     * @param m is a double and gives the mass of the particle.
     * @param lab is a String and gives the label of the particle.
     */
//Partcile's position vector, velocity vector and its label.
    public Particle3D(Vector3D X, Vector3D V, double m, String lab){
        this.x = X;
        this.v = V;
        this.mass = m;
        this.label = lab;
    }
    /**Constructs the individual vector components of position, velocity and the mass and label.
     * 
     * @param x is a double and gives the x component of the particle's position.
     * @param y is a double and gives the y component of the particle's position.
     * @param z is a double and gives the z component of the particle's position.
     * @param xv is a double and gives the x component of the particle's velocity.
     * @param yv is a double and gives the y component of the particle's velocity.
     * @param zv is a double and gives the z component of the particle's velocity.
     * @param m is a double and gives the mass of the particle.
     * @param lab is a String and gives the particle's label.
     */
    
//Particle's position and velocity components and the particle label
    public Particle3D(double x, double y, double z, double xv, double yv, double zv, double m, String lab){
        this.x = new Vector3D(x,y,z);
        this.v = new Vector3D(xv,yv,zv);
        this.mass = m;
        this.label = lab;
    }
    /////////////////////////////////////////////////////////////////////////////
    
    //Method for setting the components of particle3D; velocity, position mass and the label.
    
    /**Convenient set method to set all of the components
     * 
     * @param xx is a double and sets the x position component of the instance Particle3D.
     * @param xy is a double and sets the y position component of the instance Particle3D.
     * @param xz is a double and sets the z position component of the instance Particle3D.
     * @param vx is a double and sets the x velocity component of the instance Particle3D.
     * @param vy is a double and sets the y velocity component of the instance Particle3D.
     * @param vz is a double and sets the z velocity component of the instance Particle3D.
     */
    
    public void setParticle3D(double xx, double xy, double xz, double vx, double vy, double vz ){
        this.v.setVector3D(vx,vy,vz);
        this.x.setVector3D(xx,xy,xz);
    }
    
    /**Sets the vector components for the Particle 3D
     * 
     * @param x is 3D position vector of the particle.
     * @param v is 3D velocity vector of the particle.
     */
    
    public void setParticle3D(Vector3D x, Vector3D v){
        this.v = v;
        this.x = x;
    }
    
    /**Sets the velocity vector of the particle
     * 
     * @param v is a 3D velocity vector and sets the velocity vector
     */
    
    public void setv(Vector3D v){
        this.v = v;
    }
    
    /**Sets the position vector of the particle.
     * 
     * @param x is a 3D position vector and set the position vector
     */
    public void setx(Vector3D x){
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
    
    ///////////////////////////////////////////////////////
    //Getter method for velocity, position, mass and the label of the 3D particle
    
    /**Gets the velocity vector of Particle3D 
     * 
     * @return a velocity vector instance representing the Particle3D velocity Vector.
     */
    
    public Vector3D getv(){
        return this.v;
    }
    
    /**Gets the position vector of Particle3D 
     * 
     * @return a position vector instance representing the Particle3D position Vector.
     */
    
    public Vector3D getx(){
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
    /////////////////////////////////////////////////////////////
    
    //Prints particles position in a readable form
    
    /**Returns a String which represents the position of Particle3D.
     * 
     */
    public String toString(){
        String tempLabel = this.getLabel();
        double x1_str = this.getx().getX();
        double x2_str = this.getx().getY();
        double x3_str = this.getx().getZ();
        return "<"+tempLabel+"> <"+x1_str+"> <"+x2_str+"> <"+x3_str+">";
    }
    /////////////////////////////////////////////////////////////

    //Scans values (Label, Mass, x, y, z, xvelocity, yvelocity, zvelocity) for components from a file and prints them
    
    /**Scans the values for Particle 3D (the label, mass, position vector and velocity vector)
     * 
     * @param scan gives the values scanned from the file.
     */
    public Particle3D(Scanner scan){
        this.label = scan.next();
        this.mass = scan.nextDouble();
        this.x = new Vector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
        this.v = new Vector3D(scan.nextDouble(), scan.nextDouble(), scan.nextDouble());
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Finding the particles kinetic energy
    
   /**Gets the kinetic energy of the particles
    * 
    * @return a double representing the kinetic energy of the particle using it's mass, 
    * position vector and velocity vector.
    */
    public double getKE(){
        return 0.5*this.getMass()*(this.getv().getX()*this.getv().getX()+this.getv().getY()*this.getv().getY()
                +this.getv().getZ()*this.getv().getZ());
    }
    //updating the velocity knowing the force and timestep
    
    /**Updates the velocity using known force and timestep.
     * 
     * @param timestep is a double which updates timestep
     * @param force is a vector which updates the force
     */
    public void updateVelocity(double timestep, Vector3D force){
        this.setv(Vector3D.addVector3D(this.getv(), force.scalarMul(timestep/this.getMass())));
    }
    //updating the position knowing the timestep
    
    /**Updates the position given the timestep
     * 
     * @param timestep is a double which updates the position
     */
    public void updatePosition(double timestep){
        this.setx(Vector3D.addVector3D(this.getx(),this.getv().scalarMul(timestep)));
    }
    //update position knowing the force and timestep
    /**Updates the position given the force and the timestep
     * 
     * @param timestep is a double that updates the position
     * @param force is a vector that updates the position
     */
    public void updatePosition(double timestep, Vector3D force){
        this.setx(Vector3D.addVector3D(Vector3D.addVector3D(this.getx(),this.getv().scalarMul(timestep)),force.scalarMul(0.5*timestep*timestep/this.getMass())));
    }
    //distance between two particles   
    /**Calculates the seperation vector between two particles
     * 
     * @param a is the reference Particle3D
     * @param b is the Particle3D a distance away
     * @return the Vector3D  between the two particles pointing away from a
     */
    public Vector3D seperationAway(Particle3D a, Particle3D b){
        return Vector3D.subVector3D(a.getx(),b.getx());
    }


    //Static Methods

    /**
     * gravitational force on the particle with respect to the star
     * @param  particle Particle3D on which the force is acting 
     * @param  star   Particle3D by which the attraction is happening
     * @param  bigG   [description]
     * @return        Vector3D the force that is acting on the particle
     */
    public static Vector3D gravitationalAttraction(Particle3D particle, Particle3D star, double bigG){
        return new Vector3D();
    }


    /**
     * The total energy for the 2 particle system
     * @param  particle one Particle3D in the system
     * @param  star   another Particle3D in the system
     * @param  bigG   [description]
     * @return        a double equal to the total energyy for that pair of particles
     */
    public static double totalEnergy(Particle3D particle, Particle3D star, double bigG){
        return 1.0;   
    }


    /**
     * Outputs the velocity tangential to the orbit of the particle given it is circular about the star v = sqrt(m1 + m2/ r) m1,m2 mass of the particles and r is the radius of the orbit.
     * @param  particle the orbiting particle with the outputted velocity
     * @param  star   the particle its orbitiing around 
     * @return        a Vector3D equal to the velocity of the particle to cause a cirular orbit v = sqrt(m1 + m2/ r) m1,m2 mass of the particles and r is the radius of the orbit
     */
    //public static Vector3D tangetialVelocity(Particle3D particle, Particle3D star){
    //    return new Vector3D();        
    //}    




 

}



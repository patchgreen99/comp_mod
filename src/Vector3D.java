/** 
 * Vector3D class demonstrating setters getters of vector elements and several vector operations 
 * Class variables :  private double x_com, private double y_com, private double z_com
 *
 * @author Patrick. Green, Emily. Dowd, Craig. Young
 * @version "10/2015"
 *
 */

public class Vector3D {


    /** Object vector elements
     *
     */
    
    private double x_com;
    private double y_com;
    private double z_com;


    /*
     *  default constructor sets all elements to zero
     */

    public Vector3D () {
        // default constructor sets all elements to zero
	this.setVector3D(0,0,0);
    }

    /** Copy constructor. Constructs a new Vector3D by copying the x,y and
     *  z components of another Vector3D instance.
     *
     * @param original the Vector3D to be copied
     */
    public Vector3D (Vector3D original) {
	this.setVector3D(original.getX(), original.getY(), original.getZ() );
    }

    /** Explicit constructor. Constructs a new Vector3D from explicitly
     *  given x,y,z components.
     *
     * @param x is a double and gives the x component of the new Vector3D
     * @param y is a double and gives the y component of the new vector3D
     * @param z is a double and gives the z component of the new vector3D
     */
    public Vector3D(double x,double y, double z) {
	this.setVector3D(x,y,z);
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    /*
     * Setters and getters
     */

    /** Convenient set method to set all components.
     *
     * @param x is a double and sets the x component of the instance Vector3D
     * @param y is a double and sets the y component of the instance Vector3D
     * @param z is a double and sets the z component of the instance Vector3D
     */
    public void setVector3D(double x, double y, double z) {
	this.setX(x);
	this.setY(y);
	this.setZ(z);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    // Setters provide the internal variables

    /** Sets the x component of a Vector3D
     *
     * @param x is a double and sets x component 
     */
    public void setX(double x) { this.x_com = x; }

    /** Sets the y component of a Vector3D
     *
     * @param y is a double and sets y component 
     */
    public void setY(double y) { this.y_com = y; }

    /** Sets the z component of a Vector3D
     *
     * @param z is a double and sets z component 
     */
    public void setZ(double z) { this.z_com = z; }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////    
    // Getters provide access to private internal variables
 
    /** Gets the x component of a Vector3D.
     *
     * @return a double instance representing the Vector3D' x component.
     */
    public double getX() { return this.x_com; }

    /** Gets the y component of a Vector3D.
     *
     * @return a double instance representing the Vector3D' y component.
     */
    public double getY() { return this.y_com; }

    /** Gets the z component of a Vector3D.
     *
     * @return a double instance representing the Vector3D' z component.
     */
    public double getZ() { return this.z_com; }

   
    ////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /** Returns a String representation of the Vector3D. Methods 
     * called 'toString' are automagically called when an object
     * is printed.<br>
     * If any component is positive (e.g. for (1.0,3.0,2.0)),
     * no sign will be given for the output, if any component is negative (e.g. for (2.0,-4.0,5.0)),
     * the negative method will return a sign
     *
     * @return a string representation of the Vector3D instance
     */
    public String toString() { 
      double x_str = this.getX();
      double y_str = this.getY();
      double z_str = this.getZ();

      return "("+isNegative(x_str)+Math.abs(x_str)+","+isNegative(y_str)+Math.abs(y_str)+","+isNegative(z_str)+Math.abs(z_str)+")";
      
    }

    /**
     * Returns a sign if and only if the argument is negative
     * @param arg is a double and if negative will return a - sign
     * @result returns - sign if arg is negative
     */

    public String isNegative(double arg) {
	if ( arg < 0.0 ){ return "-";}
	else{ return "";}
    }

    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    /*
     * Instance methods
     *
     * These operate on a particular instance of an object
     */

    /** Calculates the magnitude of Vector3D A.
     * 
     * @return a double representing the Vector3D' magnitude.
     */ 
    public double magnitude() { return  Math.sqrt(this.magnitudeSquared()); }

    /** Calculates the magnitude sqaured of Vector3D A.
     *
     * @return a double representing the Vector3D' magnitude squared.
     */
    public double magnitudeSquared() { return this.getX() * this.getX() + 
                                               this.getY() * this.getY() +
	                                        this.getZ() * this.getZ();}
    

    /** Calculates the result after Vector3D A is multiplied by a scalar a.
     *
     * @param a is a double multiplied by the original components to set the new components of Vector3D
     */
    public Vector3D scalarMul(double a) {
	return new Vector3D(this.getX()*a,this.getY()*a,this.getZ()*a);	
    }

    /** Calculates the result after Vector3D A is divided by a scalar a.
     *
     * @param a is a double divided by the original components to set the new components of Vector3D
     */
    public Vector3D scalarDiv(double a) {
	return new Vector3D(this.getX()/a,this.getY()/a,this.getZ()/a);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Static Methods
    
    /** Adds two Vector3D.
     *
     * @param a the first Vector3D to be added
     * @param b the second Vector3D to be added
     * @return the Vector3D sum of a and b, a+b.
     */
    public static Vector3D addVector3D(Vector3D a, Vector3D b) { 
	return new Vector3D(a.getX() + b.getX(),
			    a.getY() + b.getY(),
			    a.getZ() + b.getZ());
    }
    
    /** Subtracts Vector3d.
     *
     * @param a the subtrahend
     * @param b the subtractor
     * @return the Vector3D difference between a and b, a-b.
     */
    public static Vector3D subVector3D(Vector3D a, Vector3D b) { 
	return new Vector3D(a.getX() - b.getX(),
			    a.getY() - b.getY(),
			    a.getZ() - b.getZ());
    }
    
    /** Derives dot product of two vectors 
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return the dot product of a and b, (a.x*b.x , a.y*b.y , a.z*b.z).
     */
    public static double dotVector3D(Vector3D a, Vector3D b) {
	return a.getX()*b.getX() + a.getY()*b.getY() + a.getZ()*b.getZ();
    }
    

    /** calculates the cross product of Vector3D a and b
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return the cross product a and b, (a.getY*b.getZ-a.getZ*b.getY , a.getZ*b.getX-a.getX*b.getZ , a.getX*b.getY-a.getY*b.getX).
     */
    public static Vector3D crossVector3D(Vector3D a, Vector3D b) {
	return new Vector3D(a.getY()*b.getZ()-a.getZ()*b.getY() , a.getZ()*b.getX()-a.getX()*b.getZ() , a.getX()*b.getY()-a.getY()*b.getX());
    }

   /** Returns the Boolean true if Vector3D a and b are equal to 4 decimal places
     *
     * @param a the first Vector3D
     * @param b the second Vector3D
     * @return true if Vector3D' are identical to 4 decimal places (each component of the vectors are equal to 4 decimal places)
     */
    public static boolean sameVector3D(Vector3D a, Vector3D b,int dp) {
	double error = 1.0/Math.pow(10,dp);

	if (Math.abs(a.getX()-b.getX())<error && Math.abs(a.getY()-b.getY())<error && Math.abs(a.getZ()-b.getZ())<error){
	    return true; } 
	else {
	    return false;}
    }

    
}


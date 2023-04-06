/**
* @file EngineMath.cpp
*
*/

#include "EngineMath.h"

//////////////////////////////////////////////////////////////////////////
// Common math functions
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// General Utility functions
//////////////////////////////////////////////////////////////////////////
#pragma region
// Are two floating point numbers equal to each other
// Floating Point Error Safe
//
// IN:		a		The first number
//			b		The second number
//
// RETURN: TRUE iff |a-b| < Tolerance
//
// NOTE:	EPSILON is tolerance
bool IsEqual(float a, float b)
{
	// NOTE: Do not modify.
	return fabs(a - b) < EPSILON;
}

// Is a floating point value equal to zero
// Floating Point Error Safe
//
// IN:		a		The number to check
//
// RETURN:	TRUE iff |a| < Tolerance
//
// NOTE:	Tolerance set by EPSILON
bool IsZero(float a)
{
	// NOTE: Do not modify
	return (fabs(a))<EPSILON;
}

// RETURN: MAX of two numbers
float Max(float a, float b)
{
	// NOTE: Do not modify.
	return (a > b) ? a : b;
}

// RETURN: MIN of two numbers
float Min(float a, float b)
{
	// NOTE: Do not modify.
	return (a < b) ? a : b;
}

// RETURN: Converts input to radian measure
float Degrees_To_Radians(float Deg)
{
	// NOTE: Do not modify.
	return Deg * PI / 180.0f;
}

// RETURN: Converts input to degree measure
float Radians_To_Degrees(float Rad)
{
	// NOTE: Do not modify.
	return Rad * 180.0f / PI;
}

float cot(float rad)
{
	return 1 / tanf(rad);
	//return cos(rad) / sin(rad);
}

float lerp(float a, float b, float ratio)
{
	return a + ratio * (b - a);
}

// CAN REMOVE THIS FUNCTION
bool is_point_on_sphere(TVECTOR point, TVECTOR center, float radius)
{
	TVECTOR distance_v = Vector_Sub(point, center);
	return (Vector_Dot(distance_v, distance_v) == radius * radius);
}

float clamp(const float& val, const float& min, const float& max)
{
	return (val - max >= EPSILON) ? max : (val - min <= EPSILON) ? min : val;
}

float saturate(const float& val)
{
	return clamp(val, 0.0f, 1.0f);
}

TVECTOR Vector_Saturate(const TVECTOR& v)
{
	TVECTOR t;
	t.x = saturate(v.x);
	t.y = saturate(v.y);
	t.z = saturate(v.z);
	t.w = saturate(v.w);

	return t;
}
#pragma endregion

////////////////////////////////////////////////////////////////////////
// Linear Algebra Functions Day 1
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Vector Functions
//////////////////////////////////////////////////////////////////////////
#pragma region VECTOR FUNCTIONS
// Check if two TVECTOR's are equal to each other
//
// IN:		v		First Vector
//			w		Second Vector
//
// RETURN:  True if v==w, False otherwise
//
// NOTE:	Use's all four components
//			Should be floating point error safe.
bool Vector_IsEqual(TVECTOR v, TVECTOR w)
{
	if (IsEqual(v.w, w.w))
	{
		if (IsEqual(v.x, w.x))
		{
			if (IsEqual(v.y, w.y))
			{
				if (IsEqual(v.z, w.z))
				{
					return true;
				}
				else
					return false;
			}
			else
				return false;
		}
		else
			return false;
	}
	else
		return false;
}

// ADD two TVECTOR's togother
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v + w
//
// NOTE:	Use's all four components
TVECTOR Vector_Add(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	temp.x = v.x + w.x;
	temp.y = v.y + w.y;
	temp.z = v.z + w.z;
	temp.w = v.w + w.w;

	return temp;
}

// SUBTRACT one TVECTOR from another
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v - w
//
// NOTE:	Use's all four components
TVECTOR Vector_Sub(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	temp.x = v.x - w.x;
	temp.y = v.y - w.y;
	temp.z = v.z - w.z;
	temp.w = v.w - w.w;

	return temp;
}

// MULTIPLY all four components of a TVECTOR by a scalar
//
// IN:		v		The vector to scale
//			s		The value to scale by
//
// RETURN:  s * v
TVECTOR Vector_Scalar_Multiply(TVECTOR v, float s)
{
	TVECTOR temp;
	temp.x = s * v.x;
	temp.y = s * v.y;
	temp.z = s * v.z;
	temp.w = s * v.w;

	return temp;
}

// NEGATE all the components of a TVECTOR
//
// IN:		v		The vector to negate
//
// RETURN:	-1 * v
//
// NOTE:	Use's all four components
TVECTOR Vector_Negate(TVECTOR v)
{
	TVECTOR temp;
	temp.x = -v.x;
	temp.y = -v.y;
	temp.z = -v.z;
	temp.w = -v.w;

	return temp;
}

// Perform a Dot Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (DOT) w
//
// NOTE:	Use's all four components
float Vector_Dot(TVECTOR v, TVECTOR w)
{
	float dot_Product = v.x * w.x + v.y * w.y + v.z * w.z + v.w * w.w;

	return dot_Product;
}

// Perform a Cross Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (CROSS) w
//
// NOTE:	The w-component of each vector is not used.
//			The resultant vector will have a w-component of zero.
TVECTOR Vector_Cross(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	temp.x = v.y * w.z - v.z * w.y;
	temp.y = -1 * (v.x * w.z - v.z * w.x);
	temp.z = v.x * w.y - v.y * w.x;
	temp.w = 0;

	return temp;
}

// Find the squared length of a TVECTOR
//
// IN:		v		The vector to find the squared length of
//
// RETURN:	Squared Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_LengthSq(TVECTOR v)
{
	float lengthSq = Vector_Dot(v, v);// v.x* v.x + v.y * v.y + v.z * v.z + v.w * v.w;
	return lengthSq;
}

// Find the length of a TVECTOR
//
// IN:		v		The vector to find the length of
//
// RETURN:	Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_Length(TVECTOR v)
{
	float length = sqrtf(Vector_LengthSq(v));
	return length;
}

// Normalize a TVECTOR
//
// IN:		v		The vector to normalize
//
// RETURN:	Normalized version of v
//
// NOTE:	Use's all four components
TVECTOR Vector_Normalize(TVECTOR v)
{
	float length = Vector_Length(v);
	if (IsZero(length))
	{
		TVECTOR zeroV;
		zeroV.x = 0;
		zeroV.y = 0;
		zeroV.z = 0;
		zeroV.w = 0;
		return zeroV;
	}
	else
	{
		TVECTOR normalizedV;
		normalizedV.x = v.x / length;
		normalizedV.y = v.y / length;
		normalizedV.z = v.z / length;
		normalizedV.w = v.w / length;
		return normalizedV;
	}
}

// Makes a TVECTOR's w-component normalized
//
// IN:		v		The vector (point object) to homogenise
//
// RETURN:	The homogenised vector (point)
//
// NOTE:	If the w-component of the vector is 0 then the
//			function will return a zero vector with a w-component
//			of 0.
TVECTOR Vector_Homogenise(TVECTOR v)
{
	if (IsZero(v.w))
	{
		TVECTOR zeroV;
		zeroV.x = 0;
		zeroV.y = 0;
		zeroV.z = 0;
		zeroV.w = 0;
		return zeroV;
	}
	else
	{
		TVECTOR h_Vector;
		h_Vector.x = v.x / v.w;
		h_Vector.y = v.y / v.w;
		h_Vector.z = v.z / v.w;
		h_Vector.w = v.w / v.w;
		return h_Vector;
	}
}

// Get a TVECTOR made from the maximun components of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A maximized vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Maximize(TVECTOR v, TVECTOR w)
{
	TVECTOR max_V;
	max_V.x = Max(v.x, w.x);
	max_V.y = Max(v.y, w.y);
	max_V.z = Max(v.z, w.z);
	max_V.w = Max(v.w, w.w);
	return max_V;
}

// Get a TVECTOR made from the minimum components of two TVECTOR's
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A minimum vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Minimize(TVECTOR v, TVECTOR w)
{
	TVECTOR min_V;
	min_V.x = Min(v.x, w.x);
	min_V.y = Min(v.y, w.y);
	min_V.z = Min(v.z, w.z);
	min_V.w = Min(v.w, w.w);
	return min_V;
}

// Get a TVECTOR made from the average of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A vector made from the average of two vectors
//
// NOTE:	Use's all four components
TVECTOR Vector_Average(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	temp.x = 0.5f * (v.x + w.x);
	temp.y = 0.5f * (v.y + w.y);
	temp.z = 0.5f * (v.z + w.z);
	temp.w = 0.5f * (v.w + w.w);
	return temp;
}

// Find the angle between two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:  The angle in degrees between the two vectors
//
// NOTE:	If either vector is a zero vector then the return
//			value will be 0.
float Vector_AngleBetween(TVECTOR v, TVECTOR w)
{
	if ((IsZero(v.x) && IsZero(v.y) && IsZero(v.z) && IsZero(v.w)) || (IsZero(w.x) && IsZero(w.y) && IsZero(w.z) && IsZero(w.w)))
	{
		return 0;
	}
	else
	{
		float rad = acosf(Vector_Dot(v, w) / (Vector_Length(v) * Vector_Length(w)));
		float angle_btw = Radians_To_Degrees(rad);
		return angle_btw;
	}
}

// Get the distance one TVECTOR points in the direction of another
// TVECTOR
//
// IN:		v		The first vector
//			w		The direction of the component
//
// RETURN:	The distance that v points in the direction of w.
//
// NOTE:	If w or v is a zero vector then the return value is zero.
float Vector_Component(TVECTOR v, TVECTOR w)
{
	if ((IsZero(w.x) && IsZero(w.y) && IsZero(w.z) && IsZero(w.w)))
	{
		return 0;
	}
	else
	{
		float comp = Vector_Dot(v, Vector_Normalize(w));
		return comp;
	}
}

// Get the TVECTOR that represents v projected on w.
//
// IN:		v		The first vector
//			w		The direction of the projection
//
// RETURN:	The projection of v onto w
//
// NOTE:	If w or v is a zero vector then the return value is zero.
TVECTOR Vector_Project(TVECTOR v, TVECTOR w)
{
	if ((IsZero(v.x) && IsZero(v.y) && IsZero(v.z) && IsZero(v.w)) || (IsZero(w.x) && IsZero(w.y) && IsZero(w.z) && IsZero(w.w)))
	{
		TVECTOR zeroV;
		zeroV.x = 0;
		zeroV.y = 0;
		zeroV.z = 0;
		zeroV.w = 0;
		return zeroV;
	}
	else
	{
		TVECTOR proj = Vector_Scalar_Multiply(Vector_Normalize(w), Vector_Component(v, w));
		return proj;
	}
}
#pragma endregion

////////////////////////////////////////////////////////////////////////
// Functions Lab  #2
///////////////////////////////////////////////////////////////////////

// Get the reflection of v across w
//
// IN:		v		The vector to reflect
//			w		The "axis" to reflect across
//
// RETURN:	v reflected across w
//
// NOTE:	If w is a zero vector then return -v.
TVECTOR Vector_Reflect(TVECTOR v, TVECTOR w)
{
	//auto reflect = [](TVECTOR v, TVECTOR w) {
	//	return Vector_Sub(Vector_Scalar_Multiply(w, 2.0f * Vector_Dot(v, w)), v);
	//};
	if (IsZero(w.x) && IsZero(w.y) && IsZero(w.z))
	{
		return Vector_Negate(v);
	}
	else if (Vector_Dot(v, w) > 0)
	{
		TVECTOR temp;
		temp = Vector_Scalar_Multiply(Vector_Project(Vector_Negate(v), w), 2);
		temp = Vector_Add(temp, v);
		return Vector_Negate(temp);
	}
	else /*if (Vector_Dot(v, w) < 0)*/
	{
		TVECTOR temp;
		temp = Vector_Scalar_Multiply(Vector_Project(v, w), 2);
		temp = Vector_Sub(temp, v);
		return temp;
	}
}

TVECTOR Vector_Lerp(TVECTOR source, TVECTOR destination, float ratio)
{
	TVECTOR v;
	v.x = lerp(source.x, destination.x, ratio);
	v.y = lerp(source.y, destination.y, ratio);
	v.z = lerp(source.z, destination.z, ratio);
	v.x = source.w;

	return v;
}

//////////////////////////////////////////////////////////////////////////
// Matrix Functions
//////////////////////////////////////////////////////////////////////////

// Get a [0] matrix
//
// RETURN: A 0 4x4 matrix
TMATRIX Matrix_Zero(void)
{
	TMATRIX zeroM;
	// col 1		// col 2		// col 3		// col 4
	zeroM._e11 = 0; zeroM._e12 = 0; zeroM._e13 = 0; zeroM._e14 = 0;	// row 1	 
	zeroM._e21 = 0; zeroM._e22 = 0; zeroM._e23 = 0; zeroM._e24 = 0; // row 2
	zeroM._e31 = 0; zeroM._e32 = 0; zeroM._e33 = 0; zeroM._e34 = 0; // row 3
	zeroM._e41 = 0; zeroM._e42 = 0; zeroM._e43 = 0; zeroM._e44 = 0; // row 4

	return zeroM;
}

// Get a [I] matrix
//
// RETURN: A 4x4 Identity matrix
TMATRIX Matrix_Identity(void)
{
	TMATRIX m = { 1, 0, 0, 0,
				  0, 1, 0, 0, 
				  0, 0, 1, 0, 
				  0, 0, 0, 1 };
	return m;
}

// Get a translation matrix
//
// IN:		x		Amount of translation in the x direction
//			y		Amount of translation in the y direction
//			z		Amount of translation in the z direction
//
// RETURN:	The translation matrix
TMATRIX Matrix_Create_Translation(float x, float y, float z)
{
	// TODO LAB 2: Replace with your implementation.
	TMATRIX m = { 1, 0, 0, x,
				  0, 1, 0, y, 
				  0, 0, 1, z,
				  0, 0, 0, 1 };
	return m;
}

// Create a scale matrix
//
// IN:		x		Amount to scale in the x direction
//			y		Amount to scale in the y direction
//			z		Amount to scale in the z direction
//
// RETURN:	The scale matrix
TMATRIX Matrix_Create_Scale(float x, float y, float z)
{
	TMATRIX scaled_M = { x, 0, 0, 0,
						 0, y, 0, 0,
						 0, 0, z, 0,
						 0, 0, 0, 1 };
	return scaled_M;
}

// Get a rotation matrix for rotation about the x-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A X-Rotation Matrix
TMATRIX Matrix_Create_Rotation_X(float Deg)
{
	float rad = Degrees_To_Radians(Deg);
	TMATRIX m = { 1, 0, 0, 0,
				  0, cos(rad), -1 * sin(rad), 0,
				  0, sin(rad), cos(rad), 0,
				  0, 0, 0, 1 };
	return m;
}

// Get a rotation matrix for rotation about the y-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Y-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Y(float Deg)
{
	float rad = Degrees_To_Radians(Deg);
	TMATRIX m = { cos(rad), 0, sin(rad), 0,
				  0, 1, 0, 0,
				  -1 * sin(rad), 0, cos(rad), 0,
				  0, 0, 0, 1 };
	return m;
}

// Get a rotation matrix for rotation about the z-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Z-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Z(float Deg)
{
	float rad = Degrees_To_Radians(Deg);
	TMATRIX m = { cos(rad), -1 * sin(rad), 0, 0,
				  sin(rad), cos(rad), 0, 0,
				  0, 0, 1, 0,
				  0, 0, 0, 1 };
	return m;
}

// ADD two matrices together
//
// IN:		m		The first matrix
//			n		The second matrix
//
// RETURN: m + n
TMATRIX Matrix_Matrix_Add(TMATRIX m, TMATRIX n)
{
	TMATRIX M;
	// col 1				  // col 2				    // col 3				  // col 4
	M._e11 = m._e11 + n._e11; M._e12 = m._e12 + n._e12; M._e13 = m._e13 + n._e13; M._e14 = m._e14 + n._e14;	// row 1	 
	M._e21 = m._e21 + n._e21; M._e22 = m._e22 + n._e22; M._e23 = m._e23 + n._e23; M._e24 = m._e24 + n._e24; // row 2
	M._e31 = m._e31 + n._e31; M._e32 = m._e32 + n._e32; M._e33 = m._e33 + n._e33; M._e34 = m._e34 + n._e34; // row 3
	M._e41 = m._e41 + n._e41; M._e42 = m._e42 + n._e42; M._e43 = m._e43 + n._e43; M._e44 = m._e44 + n._e44; // row 4

	return M;
}

// SUBTRACT two matrices
//
// IN:		m		The first matrix (left hand side)
//			n		The second matrix (right hand side)
//
// RETURN: m - n
TMATRIX Matrix_Matrix_Sub(TMATRIX m, TMATRIX n)
{
	TMATRIX M;
	// col 1				  // col 2				    // col 3				  // col 4
	M._e11 = m._e11 - n._e11; M._e12 = m._e12 - n._e12; M._e13 = m._e13 - n._e13; M._e14 = m._e14 - n._e14;	// row 1	 
	M._e21 = m._e21 - n._e21; M._e22 = m._e22 - n._e22; M._e23 = m._e23 - n._e23; M._e24 = m._e24 - n._e24; // row 2
	M._e31 = m._e31 - n._e31; M._e32 = m._e32 - n._e32; M._e33 = m._e33 - n._e33; M._e34 = m._e34 - n._e34; // row 3
	M._e41 = m._e41 - n._e41; M._e42 = m._e42 - n._e42; M._e43 = m._e43 - n._e43; M._e44 = m._e44 - n._e44; // row 4

	return M;
}

// Multiply a matrix by a scalar
//
// IN:		m		The matrix to be scaled (right hand side)
//			s		The value to scale by   (left hand side)
//
// RETURN:	The matrix formed by s*[m]
TMATRIX Matrix_Scalar_Multiply(TMATRIX m, float s)
{
	TMATRIX M;
	// col 1			 // col 2			  // col 3			   // col 4
	M._e11 = s * m._e11; M._e12 = s * m._e12; M._e13 = s * m._e13; M._e14 = s * m._e14;	// row 1	 
	M._e21 = s * m._e21; M._e22 = s * m._e22; M._e23 = s * m._e23; M._e24 = s * m._e24; // row 2
	M._e31 = s * m._e31; M._e32 = s * m._e32; M._e33 = s * m._e33; M._e34 = s * m._e34; // row 3
	M._e41 = s * m._e41; M._e42 = s * m._e42; M._e43 = s * m._e43; M._e44 = s * m._e44; // row 4
	return M;
}

// Negate a matrix
//
// IN:		m		The matrix to negate
//
// RETURN:  The negation of m
TMATRIX Matrix_Negate(TMATRIX m)
{
	TMATRIX negM;
	// col 1				 // col 2				  // col 3				   // col 4
	negM._e11 = -1 * m._e11; negM._e12 = -1 * m._e12; negM._e13 = -1 * m._e13; negM._e14 = -1 * m._e14; // row 1	 
	negM._e21 = -1 * m._e21; negM._e22 = -1 * m._e22; negM._e23 = -1 * m._e23; negM._e24 = -1 * m._e24; // row 2
	negM._e31 = -1 * m._e31; negM._e32 = -1 * m._e32; negM._e33 = -1 * m._e33; negM._e34 = -1 * m._e34; // row 3
	negM._e41 = -1 * m._e41; negM._e42 = -1 * m._e42; negM._e43 = -1 * m._e43; negM._e44 = -1 * m._e44; // row 4
	return negM;
}

// Transpose a matrix
//
// IN:		m		The matrix to transpose
//
// RETURN:	The transpose of m
TMATRIX Matrix_Transpose(TMATRIX m)
{
	TMATRIX M;
	// col 1		 // col 2		  // col 3		   // col 4
	M._e11 = m._e11; M._e12 = m._e21; M._e13 = m._e31; M._e14 = m._e41;	// row 1	 
	M._e21 = m._e12; M._e22 = m._e22; M._e23 = m._e32; M._e24 = m._e42; // row 2
	M._e31 = m._e13; M._e32 = m._e23; M._e33 = m._e33; M._e34 = m._e43; // row 3
	M._e41 = m._e14; M._e42 = m._e24; M._e43 = m._e34; M._e44 = m._e44; // row 4
	return M;
}

// Multipy a matrix and a vector
//
// IN:		m		The matrix (left hand side)
//			v		The vector (right hand side)
//
// RETURN:	[m]*v
TVECTOR Matrix_Vector_Multiply(TMATRIX m, TVECTOR v)
{
	TVECTOR V;
	// vector * row
	V.x = m._e11 * v.x + m._e12 * v.y + m._e13 * v.z + m._e14 * v.w;
	V.y = m._e21 * v.x + m._e22 * v.y + m._e23 * v.z + m._e24 * v.w;
	V.z = m._e31 * v.x + m._e32 * v.y + m._e33 * v.z + m._e34 * v.w;
	V.w = m._e41 * v.x + m._e42 * v.y + m._e43 * v.z + m._e44 * v.w;
	return V;
}

// Multipy a vector and a matrix
//
// IN:		v		The vector ( left hand side)
//			m		The matrix (right hand side)
//
// RETURN:	v*[m]
TVECTOR Vector_Matrix_Multiply(TVECTOR v, TMATRIX m)
{
	TVECTOR V;
	// vector * column
	V.x = v.x * m._e11 + v.y * m._e21 + v.z * m._e31 + v.w * m._e41;
	V.y = v.x * m._e12 + v.y * m._e22 + v.z * m._e32 + v.w * m._e42;
	V.z = v.x * m._e13 + v.y * m._e23 + v.z * m._e33 + v.w * m._e43;
	V.w = v.x * m._e14 + v.y * m._e24 + v.z * m._e34 + v.w * m._e44;
	return V;
}
// Multiply a matrix by a matrix
//
// IN:		m		First Matrix (left hand side)
//			n		Second Matrix (right hand side)
//
// RETURN:	[m]*[n]
TMATRIX Matrix_Matrix_Multiply(TMATRIX m, TMATRIX n)
{
	TMATRIX M;
	// row 1
	M._e11 = m._e11 * n._e11 + m._e12 * n._e21 + m._e13 * n._e31 + m._e14 * n._e41; // m_row 1 * n_col 1
	M._e12 = m._e11 * n._e12 + m._e12 * n._e22 + m._e13 * n._e32 + m._e14 * n._e42; // m_row 1 * n_col 2
	M._e13 = m._e11 * n._e13 + m._e12 * n._e23 + m._e13 * n._e33 + m._e14 * n._e43; // m_row 1 * n_col 3
	M._e14 = m._e11 * n._e14 + m._e12 * n._e24 + m._e13 * n._e34 + m._e14 * n._e44; // m_row 1 * n_col 4
	// row 2
	M._e21 = m._e21 * n._e11 + m._e22 * n._e21 + m._e23 * n._e31 + m._e24 * n._e41;	// m_row 2 * n_col 1
	M._e22 = m._e21 * n._e12 + m._e22 * n._e22 + m._e23 * n._e32 + m._e24 * n._e42;	// m_row 2 * n_col 2
	M._e23 = m._e21 * n._e13 + m._e22 * n._e23 + m._e23 * n._e33 + m._e24 * n._e43;	// m_row 2 * n_col 3
	M._e24 = m._e21 * n._e14 + m._e22 * n._e24 + m._e23 * n._e34 + m._e24 * n._e44; // m_row 2 * n_col 4
	// row 3
	M._e31 = m._e31 * n._e11 + m._e32 * n._e21 + m._e33 * n._e31 + m._e34 * n._e41;	// m_row 3 * n_col 1
	M._e32 = m._e31 * n._e12 + m._e32 * n._e22 + m._e33 * n._e32 + m._e34 * n._e42; // m_row 3 * n_col 2
	M._e33 = m._e31 * n._e13 + m._e32 * n._e23 + m._e33 * n._e33 + m._e34 * n._e43; // m_row 3 * n_col 3
	M._e34 = m._e31 * n._e14 + m._e32 * n._e24 + m._e33 * n._e34 + m._e34 * n._e44; // m_row 3 * n_col 4
	// row 4
	M._e41 = m._e41 * n._e11 + m._e42 * n._e21 + m._e43 * n._e31 + m._e44 * n._e41; // m_row 4 * n_col 1
	M._e42 = m._e41 * n._e12 + m._e42 * n._e22 + m._e43 * n._e32 + m._e44 * n._e42; // m_row 4 * n_col 2
	M._e43 = m._e41 * n._e13 + m._e42 * n._e23 + m._e43 * n._e33 + m._e44 * n._e43; // m_row 4 * n_col 3
	M._e44 = m._e41 * n._e14 + m._e42 * n._e24 + m._e43 * n._e34 + m._e44 * n._e44; // m_row 4 * n_col 4
	return M;
}

////////////////////////////////////////////////////////////////////////
// Matrix Functions Lab # 3
///////////////////////////////////////////////////////////////////////

// HELPER FUNCTION  *** NOT GRADED, ONLY SUGGESTED ***
// USE THIS FUNCTION TO FIND THE DETERMINANT OF A 3*3
// MATRIX. IT CAN BE USED IN THE MATRIX DETERMINANT
// AND MATRIX INVERSE FUNCTIONS BELOW
// 
// RETURN:	The determinant of a 3x3 matrix
float Matrix_Determinant(float e_11,float e_12,float e_13,
						 float e_21,float e_22,float e_23,
						 float e_31,float e_32,float e_33)
{
	float det_3x3 = e_11 * (e_22 * e_33 - e_32 * e_23) - e_12 * (e_21 * e_33 - e_31 * e_23) + e_13 * (e_21 * e_32 - e_31 * e_22);
	return det_3x3;
}

// Get the determinant of a matrix
//
// IN:		m		The ONE!
//
// RETURN:	It's determinant
float Matrix_Determinant(TMATRIX m)
{
	float det_11 = m._e11 * Matrix_Determinant(m._e22, m._e23, m._e24, m._e32, m._e33, m._e34, m._e42, m._e43, m._e44);
	float det_12 = m._e12 * Matrix_Determinant(m._e21, m._e23, m._e24, m._e31, m._e33, m._e34, m._e41, m._e43, m._e44);
	float det_13 = m._e13 * Matrix_Determinant(m._e21, m._e22, m._e24, m._e31, m._e32, m._e34, m._e41, m._e42, m._e44);
	float det_14 = m._e14 * Matrix_Determinant(m._e21, m._e22, m._e23, m._e31, m._e32, m._e33, m._e41, m._e42, m._e43);
	float Det = det_11 - det_12 + det_13 - det_14;
	return Det;
}

// Get the inverse of a matrix
//
// IN:		m		The matrix to inverse
//
// RETURN:	The Inverse of [m]
//
// NOTE: Returns the matrix itself if m is not invertable.
TMATRIX Matrix_Inverse(TMATRIX m)
{
	if (Matrix_Determinant(m) != 0)
	{
		TMATRIX Adj_M;
		Adj_M._e11 = Matrix_Determinant(m._e22, m._e23, m._e24, m._e32, m._e33, m._e34, m._e42, m._e43, m._e44);
		Adj_M._e12 = -1 * Matrix_Determinant(m._e21, m._e23, m._e24, m._e31, m._e33, m._e34, m._e41, m._e43, m._e44);
		Adj_M._e13 = Matrix_Determinant(m._e21, m._e22, m._e24, m._e31, m._e32, m._e34, m._e41, m._e42, m._e44);
		Adj_M._e14 = -1 * Matrix_Determinant(m._e21, m._e22, m._e23, m._e31, m._e32, m._e33, m._e41, m._e42, m._e43);

		Adj_M._e21 = -1 * Matrix_Determinant(m._e12, m._e13, m._e14, m._e32, m._e33, m._e34, m._e42, m._e43, m._e44);
		Adj_M._e22 = Matrix_Determinant(m._e11, m._e13, m._e14, m._e31, m._e33, m._e34, m._e41, m._e43, m._e44);
		Adj_M._e23 = -1 * Matrix_Determinant(m._e11, m._e12, m._e14, m._e31, m._e32, m._e34, m._e41, m._e42, m._e44);
		Adj_M._e24 = Matrix_Determinant(m._e11, m._e12, m._e13, m._e31, m._e32, m._e33, m._e41, m._e42, m._e43);

		Adj_M._e31 = Matrix_Determinant(m._e12, m._e13, m._e14, m._e22, m._e23, m._e24, m._e42, m._e43, m._e44);
		Adj_M._e32 = -1 * Matrix_Determinant(m._e11, m._e13, m._e14, m._e21, m._e23, m._e24, m._e41, m._e43, m._e44);
		Adj_M._e33 = Matrix_Determinant(m._e11, m._e12, m._e14, m._e21, m._e22, m._e24, m._e41, m._e42, m._e44);
		Adj_M._e34 = -1 * Matrix_Determinant(m._e11, m._e12, m._e13, m._e21, m._e22, m._e23, m._e41, m._e42, m._e43);

		Adj_M._e41 = -1 * Matrix_Determinant(m._e12, m._e13, m._e14, m._e22, m._e23, m._e24, m._e32, m._e33, m._e34);
		Adj_M._e42 = Matrix_Determinant(m._e11, m._e13, m._e14, m._e21, m._e23, m._e24, m._e31, m._e33, m._e34);
		Adj_M._e43 = -1 * Matrix_Determinant(m._e11, m._e12, m._e14, m._e21, m._e22, m._e24, m._e31, m._e32, m._e34);
		Adj_M._e44 = Matrix_Determinant(m._e11, m._e12, m._e13, m._e21, m._e22, m._e23, m._e31, m._e32, m._e33);

		TMATRIX inv_M = Matrix_Transpose(Matrix_Scalar_Multiply(Adj_M, 1 / Matrix_Determinant(m)));
		return inv_M;
	}
	else
		return  m;
}

// Get a projection matrix for prepping perspective
//
// IN:		FOV		
//			nearPlane
//			farPlane
//
// RETURN:	A Projection Matrix
TMATRIX Matrix_Create_Projection(unsigned int FOV, float nearPlane, float farPlane, int width, int height)
{
	float Yscale = cot(Degrees_To_Radians(static_cast<float>(FOV >> 1)));
	float Xscale = Yscale * (static_cast<float>(height) / width);
	TMATRIX ProjectionMatrix = { Xscale, 0, 0, 0,
							0, Yscale, 0, 0,
							0, 0, farPlane / (farPlane - nearPlane), 1,
							0, 0, -(farPlane * nearPlane) / (farPlane - nearPlane), 0 };
	return ProjectionMatrix;
}
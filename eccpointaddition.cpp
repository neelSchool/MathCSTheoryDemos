/*
you have a curve E in Friars Cross form y^2 = x^3 + Ax + B
you have points P (x1,y1) and Q (x2,y2) on E which need to be added.
you need to return the answer which is R'.
The process is that you draw a line through point P and Q which will then intersect E at a 3rd point which is R.
The two edge cases can be done later but this is the general case.
*/
#include <iostream>
using namespace std;
void pointadd(float px, float py, float qx, float qy, float A, float B) {
    //calculate line L of form y = mx + c
    float m = (qy - py) / (qx - px);
    float c = py - m * px;
    // sub y = mx+c into y^2=x^3+Ax+B
    // (mx+c)^2 = m^2x^2 + c^2 + 2mx = x^3 + Ax + B
    // x^3 - m^2x^2 + (A-2m)x + (B-c^2) = 0
    // the x coords of p and q give 2 out of 3 roots
    // so x coory of R will be the third root.
    // using vietas formula r1 + r2 + r3 = -b/a.
    // since a (coeff of x^3) is 1 and b is -m^2.
    // r3 = m^2 - px - qx
    float rx = (m*m) - px - qx;
    float rprimex = rx;
    float ry = m*rx + c;
    float rprimey = -1 * ry;
    cout << "(" << rprimex << "," << rprimey << ")"; 
}
int main() {
    // example settings
    // Curve E is y^2 = x^3 + 0x + 1, Point P is (-1,0) and Q is (0,1)
    // so R is (2,3) and R' is (2,-3)
    // lets see if i can build an algorithm to get this to work.
    float px = -1, py = 0, qx = 0, qy = 1; // p and q points
    float A = 0, B = 1; //curve coefficients
    pointadd(px,py,qx,qy,A,B);
}
// demo by Neel Mehendale 
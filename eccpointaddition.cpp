/*
you have a curve E in Friars Cross form y^2 = x^3 + Ax + B
you have points P (x1,y1) and Q (x2,y2) on E which need to be added.
you need to return the answer which is R'.
The process is that you draw a line through point P and Q which will then intersect E at a 3rd point which is R.
*/
// Note: All examples given were verified using DESMOS Graphing Calculator and Casio FX-CG50. They are all whole number examples to serve as an easy to understand proof of concept.
// An implementation where this is applied properly to Finite Fields will be done later as a separate file.
#include <iostream>
using namespace std;
void horizontalcase(float px, float py) {
    cout << "(" << px << "," << py << ")";
}
void verticalcase(float px, float py) {
    // If P and Q lie on a vertical line then the line hits Infinity which acts as an Identity (or zero) so P+Q=P+0=P
    // Example would be if P is (0,-1) and Q is (0,1) then P+Q = (0,-1)
    cout <<  "(" << px << "," << py << ")";
}
void tangentialcase(float px, float py, float A) {
    // If P and Q are the same point then P + Q = P + P = R' by drawing a tangent at P on E and find second intersection and reflecting on x-axis
    // Example would be if P and Q are both input as (2,3) then tangent intersects at R which is (0,-1) so R' is (0,1)
    // Find gradient of y = mx+c tangent line by differentaiting and subbing in x and y coord of P
    // dy/dx = (3x^2 / 2y) + A 

    // However if you take P as (0,-1) you get a stationary point with a horizontal case (dy/dx = 0), which is another infinity or zero like the vertical case. So P+Q=P+P=P+0=P
    float m = ((3 * px * px) / (2 * py)) + A;
    if (m == 0) {
        horizontalcase(px,py);
        return;
    }
    float c = py - (m * px);
    // now use general solution of rx = m^2 - 2px by factorising cubic
    float rx = (m * m) - (2 * px);
    float ry = (2 * rx) + c;
    float rprimex = rx, rprimey = -1 * ry;
    cout << "(" << rprimex << "," << rprimey << ")"; 
}
void pointadd(float px, float py, float qx, float qy) {
    //calculate line L of form y = mx + c
    float m = (qy - py) / (qx - px);
    float c = py - m * px;
    // sub y = mx+c into y^2=x^3+Ax+B
    // (mx+c)^2 = m^2x^2 + c^2 + 2mx = x^3 + Ax + B
    // x^3 - m^2x^2 + (A-2m)x + (B-c^2) = 0
    // the x coords of p and q give 2 out of 3 roots
    // so x coord of R will be the third root.
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
    /*// example settings
    // Curve E is y^2 = x^3 + 0x + 1, Point P is (-1,0) and Q is (0,1)
    // so R is (2,3) and R' is (2,-3)
    // lets see if i can build an algorithm to get this to work.
    float px = -1, py = 0, qx = 0, qy = 1; // p and q points
    float A = 0, B = 1; //curve coefficients
    pointadd(px,py,qx,qy); */
    // B is not actually required in any of these cases but we take input for clarity and in case there are more edge cases which I am not yet aware of that do need B.
    cout << "Enter space separated px, py, qx, qy, A, B." << "\n";
    float px, py, qx, qy, A, B;
    cin >> px >> py >> qx >> qy >> A >> B;
    if (px == qx && py == qy) {
        tangentialcase(px,py, A); // can lead to standard tangent case or stationary horizontal case
    } else if (px == qx) {
        verticalcase(px,py);
    } else {
        pointadd(px,py,qx,qy); // general case
    }
}
// demo by Neel Mehendale 

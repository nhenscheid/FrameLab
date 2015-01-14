#include<iostream>
#include"find_l_3d.h"


int main()
{
    float x1,x2,y1,y2,z1,z2,dx,dy,dz,x,y,z,l;
    
    x1 = 0.0;
    y1 = 0.0;
    z1 = 0.0;
    x2 = 0.0;
    y2 = 0.0;
    z2 = 1.0;
    dx = 0.25;
    dy = 0.5;
    dz = 0.25;
    x = 5.0;
    y = 5.0;
    z = 5.0;
    

    l = find_l_3d(x1,y1,z1,x2,y2,z2,dx,dy,dz,x,y,z);
    
    std::cout<<"l = "<<l<<std::endl;


    return 0;
}

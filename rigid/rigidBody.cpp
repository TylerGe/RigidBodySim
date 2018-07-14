#include "rigidBody.h"

RigidBody::RigidBody()
{
	Vector3d center(0, 20, 0);
	vertexPos[0] = Vector3d(-halfSize, -halfSize, halfSize) + center;
	vertexPos[1] = Vector3d(halfSize, -halfSize, halfSize) + center;
	vertexPos[2] = Vector3d(halfSize, -halfSize, -halfSize) + center;
	vertexPos[3] = Vector3d(-halfSize, -halfSize, -halfSize) + center;
	vertexPos[4] = Vector3d(-halfSize, halfSize, halfSize) + center;
	vertexPos[5] = Vector3d(halfSize, halfSize, halfSize) + center;
	vertexPos[6] = Vector3d(halfSize, halfSize, -halfSize) + center;
	vertexPos[7] = Vector3d(-halfSize, halfSize, -halfSize) + center;

	for (int i = 0; i < number; i++) {
		vertexPosNew[i] = vertexPos[i];
	}

	//Matrix3x3 m(0, 0, 0, 0, 0, 0, 0, 0, 0);
	Matrix3x3 m(1, 0, 0, 0, 1, 0, 0, 0, 1);
	rigidState.xposition = center;
	rigidState.quater = Quaternion(m);
	rigidState.pfmom = Vector3d(0, 0, 0);
	rigidState.lamom = Vector3d(0, 0, 0);

	Vector3d w(0, 0, 0);
	rigidStateDot.velocity = Vector3d(0, 0, 0);
	rigidStateDot.quaterA = 0.5*w*(rigidState.quater);
	rigidStateDot.force = Vector3d(0, 0, 0);
	rigidStateDot.torque = Vector3d(0, 0, 0);

}

RigidBody::~RigidBody()
{
}

void RigidBody::startFall()
{
	resetSign = false;

	bodyForce = Vector3d(0, -10, 0);
}

Matrix3x3 RigidBody::getStateRotation()
{
	Matrix3x3 m = rigidState.quater.rotation();
	return m;
}
 
Vector3d RigidBody::getStatePosition()
{
	return rigidState.xposition;
}

void RigidBody::updateFall()
{
    double hStep = 0.01;
    if(resetSign){
        rigidStateDot.velocity=0;
        resetSign=false;
    }
    statesNumInt(rigidState, rigidStateDot, rigidStateNew, hStep);
    collisionDetect(rigidState, rigidStateDot, rigidStateNew, hStep);
    
    
    rigidState = rigidStateNew;
    center = rigidState.xposition;
    for (int i=0; i<number; i++) {
        vertexPos[i] = vertexPosNew[i];
    }

}

StateDot RigidBody::F(Rigidstate& rigidState)// compute force
{
    StateDot rigidStateDot;
    
    rigidStateDot.velocity = (1.0/mass)*rigidState.pfmom;
    Matrix3x3 R = rigidState.quater.rotation();
    Matrix3x3 Iinverse = R * (Io.inv()) * (R.transpose());
    Vector3d w = Iinverse * rigidState.lamom;
    rigidStateDot.quaterA = 0.5 * w * rigidState.quater;
    
    Vector3d totalforce(0,0,0);
    for (int i=0; i<number; i++) {
        totalforce = totalforce + vertexForce[i];
    }
    totalforce = totalforce + bodyForce;
    rigidStateDot.force = totalforce ;
    
    
    Vector3d totalTorque(0,0,0);
    Vector3d torque(0,0,0);
    for (int i=0; i<number; i++) {
        torque = (vertexPos[i] - center) % vertexForce[i];
        totalTorque = totalTorque + torque;
    }
    rigidStateDot.torque = totalTorque;
    
    return rigidStateDot;
}

void RigidBody::statesNumInt(Rigidstate& rigidState, StateDot& rigidStateDot, Rigidstate& rigidStateNew, double h)
{

    StateDot K1 = F(rigidState);
    Rigidstate tmp1 = rigidState + K1*h;

    StateDot K2 = F(tmp1);
    rigidStateDot = (K1 + K2)*(0.5);
    rigidStateNew = rigidState + rigidStateDot*(h);
    
    Vector3d centerNew = rigidStateNew.xposition;
    Matrix3x3 R = rigidStateNew.quater.rotation();
    
    vertexPosNew[0] = R*(Vector3d(-halfSize,-halfSize,halfSize))+centerNew;
    vertexPosNew[1] = R*(Vector3d(halfSize,-halfSize,halfSize))+centerNew;
    vertexPosNew[2] = R*(Vector3d(halfSize,-halfSize,-halfSize))+centerNew;
    vertexPosNew[3] = R*(Vector3d(-halfSize,-halfSize,-halfSize))+centerNew;
    vertexPosNew[4] = R*(Vector3d(-halfSize,halfSize,halfSize))+centerNew;
    vertexPosNew[5] = R*(Vector3d(halfSize,halfSize,halfSize))+centerNew;
    vertexPosNew[6] = R*(Vector3d(halfSize,halfSize,-halfSize))+centerNew;
    vertexPosNew[7] = R*(Vector3d(-halfSize,halfSize,-halfSize))+centerNew;

}

void RigidBody::collisionDetect(Rigidstate& rigidState, StateDot& rigidStateDot, Rigidstate& rigidStateNew, double h)
{
    Vector3d floorNormal={0,1,0};
    
    
    for (int i=0; i<number; i++) {
        double axbyczd = vertexPos[i]*floorNormal;
        if(axbyczd<DepthEpsilon){
            h=h*0.5;
            statesNumInt(rigidState, rigidStateDot, rigidStateNew, h);
            rigidState = rigidStateNew;
            center = rigidState.xposition;
            for (int i=0; i<number; i++) {
                vertexPos[i] = vertexPosNew[i];
            }
        }
        if(axbyczd>=thr){
            statesNumInt(rigidState, rigidStateDot, rigidStateNew, h);
        }
        if(axbyczd<thr){
            unsigned int Counter=0;
            if(Counter<5)
            {
                ResolveCollisions(i);
                Counter++;
            }
            else{
                resetSign=true;
                return;
            }
            statesNumInt(rigidState, rigidStateDot, rigidStateNew, h);
        }
    }
	
}

void RigidBody::ResolveCollisions(int cornerIndex)
{
    Matrix3x3 R = rigidState.quater.rotation();
    Matrix3x3 Iinverse = R * (Io.inv()) * (R.transpose());
    Vector3d w = Iinverse * rigidState.lamom;
    Vector3d r = vertexPos[cornerIndex] - center;

    double Vn =  Vpn * (rigidStateDot.velocity + w%r);
    double j = -1*(1+Cr)*Vn / ( 1.0/mass + Vpn* (Iinverse*(r%Vpn)%r) );
    Vector3d J = j * Vpn;

    rigidState.pfmom = rigidState.pfmom + J ;
    rigidState.lamom = rigidState.lamom + r%J ;
}

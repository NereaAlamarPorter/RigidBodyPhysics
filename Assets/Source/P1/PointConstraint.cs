using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;

/// <summary>
/// Basic point constraint between two rigid bodies.
/// </summary>
public class PointConstraint : MonoBehaviour, IConstraint
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public PointConstraint()
    {
        Manager = null;
    }

    #region EditorVariables

    public float Stiffness;

    public RigidBody bodyA;
    public RigidBody bodyB;

    #endregion

    #region OtherVariables

    int index;
    private PhysicsManager Manager;

    protected Vector3 pointA;
    protected Vector3 pointB;

    #endregion

    #region MonoBehaviour

    // Update is called once per frame
    void Update()
    {
        // Compute the average position
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
        Vector3 pos = 0.5f * (posA + posB);

        // Apply the position
        Transform xform = GetComponent<Transform>();
        xform.position = pos;
    }

    #endregion

    #region IConstraint

    public void Initialize(int ind, PhysicsManager m)
    {
        index = ind;
        Manager = m;

        // Initialize local positions. We assume that the object is connected to a Sphere mesh.
        Transform xform = GetComponent<Transform>();
        if (xform == null)
        {
            System.Console.WriteLine("[ERROR] Couldn't find any transform to the constraint");
        }
        else
        {
            System.Console.WriteLine("[TRACE] Succesfully found transform connected to the constraint");
        }

        // Initialize kinematics
        Vector3 pos = xform.position;

        // Local positions on objects
        pointA = (bodyA != null) ? bodyA.PointGlobalToLocal(pos) : pos;
        pointB = (bodyB != null) ? bodyB.PointGlobalToLocal(pos) : pos;

    }

    public int GetNumConstraints()
    {
        return 3;
    }

    public void GetConstraints(VectorXD c)
    {
        // TO BE COMPLETED
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
        float cx, cy, cz;
        cx = (posA[0] - posB[0]);
        cy = (posA[1] - posB[1]);
        cz = (posA[2] - posB[2]);


        c.SetSubVector(index, 3, Utils.ToVectorXD(new Vector3(cx,cy,cz)));


    }

    public void GetConstraintJacobian(MatrixXD dcdx) //Gradiente de C
    {
        // TO BE COMPLETED
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;
        MatrixXD dcdPosa = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdPosb = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdQa = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdQb = MatrixXD.Build.Dense(3, 3);
        MatrixXD I = DenseMatrixXD.CreateIdentity(3);

        float cx, cy, cz;
        cx = (posA[0] - posB[0]);
        cy = (posA[1] - posB[1]);
        cz = (posA[2] - posB[2]);


        if (bodyA != null)
        {
            dcdPosa = I;
            dcdQa = -Utils.Skew(posA - bodyA.m_pos);
            for (int i = 0; i < 3; i++) {
                dcdx.SetSubMatrix(index, bodyA.index, dcdx.SubMatrix(index, 3, bodyA.index,3) +dcdPosa);
                dcdx.SetSubMatrix(index, bodyA.index+3, dcdx.SubMatrix(index, 3, bodyA.index+3, 3) + dcdQa);
            }
        }
        if (bodyB != null)
        {
            dcdPosb = -I;
            dcdQb = Utils.Skew(posB - bodyB.m_pos);
            for (int i = 0; i < 3; i++)
            {
                dcdx.SetSubMatrix(index, bodyB.index, dcdx.SubMatrix(index, 3, bodyB.index, 3) + dcdPosb);
                dcdx.SetSubMatrix(index, bodyB.index + 3, dcdx.SubMatrix(index, 3, bodyB.index + 3, 3) + dcdQb);
            }
        }
        //Esto es gracias a mi explicación super a sucio que voy a pasar a limpio a que sí Nerea? :D


    }

    public void GetForce(VectorXD force)
    {
        // TO BE COMPLETED
        Vector3 posA = (bodyA != null) ? bodyA.PointLocalToGlobal(pointA) : pointA;
        Vector3 posB = (bodyB != null) ? bodyB.PointLocalToGlobal(pointB) : pointB;

        float cx, cy, cz;
        cx = (posA[0] - posB[0]);
        cy = (posA[1] - posB[1]);
        cz = (posA[2] - posB[2]);

        MatrixXD dcdPosa = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdPosb = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdQa = MatrixXD.Build.Dense(3, 3);
        MatrixXD dcdQb = MatrixXD.Build.Dense(3, 3);
        MatrixXD I = DenseMatrixXD.CreateIdentity(3);

        if (bodyA != null)
        {
            dcdPosa = I;
            dcdQa = -Utils.Skew(posA - bodyA.m_pos);
            VectorXD faPos = -Stiffness * dcdPosa.Transpose() * Utils.ToVectorXD(new Vector3(cx, cy, cz));
            VectorXD faQ = -Stiffness * dcdQa.Transpose() * Utils.ToVectorXD(new Vector3(cx, cy, cz));
            force.SetSubVector(bodyA.index, 3, force.SubVector(bodyA.index, 3) + faPos);
            force.SetSubVector(bodyA.index + 3, 3, force.SubVector(bodyA.index + 3, 3) + faQ);
        }

        if (bodyB != null)
        {
            dcdPosb = -I;
            dcdQb = Utils.Skew(posB - bodyB.m_pos);
            VectorXD fbPos = -Stiffness * dcdPosb.Transpose() * Utils.ToVectorXD(new Vector3(cx, cy, cz));
            VectorXD fbQ = -Stiffness * dcdQb.Transpose() * Utils.ToVectorXD(new Vector3(cx, cy, cz));
            force.SetSubVector(bodyB.index, 3, force.SubVector(bodyB.index, 3) + fbPos);
            force.SetSubVector(bodyB.index + 3, 3, force.SubVector(bodyB.index + 3, 3) + fbQ);
        }
        
    }

    public void GetForceJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        // TO BE COMPLETED
    }

    #endregion

}

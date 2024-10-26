using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class XpbdFluidParticel
{
    public Vector3 x; // position
    public Vector3 p; // predicted position
    public Vector3 v; // velocity
    public float lambda;

    public XpbdFluidParticel(Vector3 x)
    {
        this.x = x;
        this.v = Vector3.zero;
    }
}

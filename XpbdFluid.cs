using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.VisualScripting;
using UnityEngine;
using static UnityEngine.ParticleSystem;




/*
     __   _______  ____  _____     ______ _       _     _ 
     \ \ / /  __ \|  _ \|  __ \   |  ____| |     (_)   | |
      \ V /| |__) | |_) | |  | |  | |__  | |_   _ _  __| |
       > < |  ___/|  _ <| |  | |  |  __| | | | | | |/ _` |
      / . \| |    | |_) | |__| |  | |    | | |_| | | (_| |
     /_/ \_\_|    |____/|_____/   |_|    |_|\__,_|_|\__,_|

 */




public class XpbdFluid : MonoBehaviour
{
    [Header("Simulation Domain")]
    public float sizeX;
    public float sizeY;


    [Header("Simulation Parameters")]
    public int numberOfParticles;
    public Vector3 gravity;
    public float smoothingRadius;
    public float restDensity;
    public float epsilon;

    public int substeps;
    private float dtf = 0.02f;
    private float dts;
    private float dts2;


    [Header("Simulation Particles")]
    public List<XpbdFluidParticel> particles = new List<XpbdFluidParticel>();
    public Mesh particleMesh;
    public Material particleMaterial;
    private List<XpbdFluidParticel>[,] neighbours;
    private int gridSizeX, gridSizeY;


    [Header("Mouse Force")]
    Vector3 mousePos;
    float strength = 0.0f;
    public float mouseRadius;
    public float mouseStrength;

    private void Start()
    {
        dts = dtf / substeps;
        dts2 = dts * dts;

        initParticles();
    }

     
    //---------------------------------------------------------------------------------------------//
    //---------------------------- Main Simulation Loop--------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    private void FixedUpdate()
    {

        for (int i = 0; i < substeps; i++)
        {
            //Predict Positions
            for (int j = 0; j < particles.Count; j++){
                Vector3 force = InteractionForce(mousePos, mouseRadius, strength, j);
                particles[j].p = particles[j].x + dts * particles[j].v + dts2 * (gravity + force);
            }

            //Initialize Neighbour Search
            initNeighbourSearch();

            //delta Lamda
            Parallel.For(0, numberOfParticles, j =>
            {
                List<XpbdFluidParticel> neighbors = getNeighbours(particles[j].p);
                float numerator = calcConstraint(j, neighbors);

                float denominator = 0;
                for(int k = 0; k < neighbors.Count; k++)
                {
                    Vector3 grad = calcGradConstraint(j, particles.IndexOf(neighbors[k]), neighbors);
                    denominator += grad.sqrMagnitude;
                }

                denominator += epsilon;
                particles[j].lambda = -numerator / denominator;
            
            });


            //delta X       
            Parallel.For(0, numberOfParticles, j =>
            {
                List<XpbdFluidParticel> neighbors = getNeighbours(particles[j].p);
                Vector3 sum = Vector3.zero;

                for (int k = 0; k < neighbors.Count; k++)
                {
                    float coeff = particles[j].lambda + neighbors[k].lambda;
                    sum += coeff * calcGradSpikyKernel(particles[j].p - neighbors[k].p, smoothingRadius);
                }

                particles[j].p += sum / restDensity;
            });
            
            solveCollisions();
            

            //Update Positions
            for (int j = 0; j < particles.Count; j++){
                particles[j].v = (particles[j].p - particles[j].x) / dts;
                particles[j].x = particles[j].p;
            }
        }
    }

    void Update()
    {
        mouseInput();
        renderParicles();
    }


    //---------------------------------------------------------------------------------------------//
    //----------------- Initialize Particles at random positions-----------------------------------//
    //---------------------------------------------------------------------------------------------//
    private void initParticles()
    {
        for (int i = 0; i < numberOfParticles; i++){
            Vector3 position = new Vector3(Random.Range(-sizeX / 2f, sizeX / 2f), Random.Range(-sizeY / 2f, sizeY / 2f), 0f);
            particles.Add(new XpbdFluidParticel(position));
        }
    }


    //---------------------------------------------------------------------------------------------//
    //------------------------- Solver Simple Collisions ------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    private void solveCollisions()
    {
        float bottom = -sizeY/2;
        float top = sizeY/2;
        float left = -sizeX/2;
        float right = sizeX/2;

        for (int j = 0; j < particles.Count; j++){
            if (particles[j].p.y < bottom)
                particles[j].p.y = bottom;

            if (particles[j].p.y > top)
                particles[j].p.y = top;

            if (particles[j].p.x < left)
                particles[j].p.x = left;

            if (particles[j].p.x > right)
                particles[j].p.x = right;         
        }
    }


    //---------------------------------------------------------------------------------------------//
    //----------------------------- Fluid Equations -----------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    private float calcConstraint(int i, List<XpbdFluidParticel> neighbors)
    {
        float density = calcDensity(i, neighbors);
        return (density / restDensity) - 1;
    }
    private float calcDensity(int i, List<XpbdFluidParticel> neighbors)
    {
        float density = 0;

        for(int j = 0; j < neighbors.Count; j++){
            density += calcPoly6Kernel(particles[i].p - neighbors[j].p, smoothingRadius);
        }
        
        return density;
    }
    private Vector3 calcGradConstraint(int i, int j, List<XpbdFluidParticel> neighbors)
    {
        if(i == j){
            Vector3 sum = Vector3.zero;
            for(int k = 0; k < neighbors.Count; k++){
                sum += calcGradSpikyKernel(particles[i].p - neighbors[k].p, smoothingRadius);
            }
       
            return sum / restDensity;
        }

        else{
            return -calcGradSpikyKernel(particles[i].p - particles[j].p, smoothingRadius) / restDensity;
        }

    }

    //---------------------------------- Kernels ------------------------------------------//
    private float calcPoly6Kernel(Vector3 dst, float h)
    {
        float dstSquared = dst.sqrMagnitude;
        float hSquared = h * h;

        if(dstSquared > hSquared || dst == Vector3.zero)
            return 0;

        float coeff = 315f / (64f * Mathf.PI * Mathf.Pow(h, 9));
        float value = Mathf.Pow(hSquared - dstSquared, 3);

        return coeff * value;
    }
    private Vector3 calcGradSpikyKernel(Vector3 dst, float h)
    {
        float dstMagnitude = dst.magnitude;

        if (dstMagnitude > h || dst == Vector3.zero)
            return Vector3.zero;

        float coeff = 45f / (Mathf.PI * Mathf.Pow(h, 6));
        float value = Mathf.Pow(h - dstMagnitude, 2);
        
        Vector3 dstNormalized = dst.normalized;
        
        return -coeff * value * dstNormalized;
    }


    //---------------------------------------------------------------------------------------------//
    //---------------------------- Neighbour Search -----------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    private void initNeighbourSearch()
    {
        gridSizeX = Mathf.CeilToInt(sizeX / smoothingRadius);
        gridSizeY = Mathf.CeilToInt(sizeY / smoothingRadius);

        neighbours = new List<XpbdFluidParticel>[gridSizeX, gridSizeY];

        for (int i = 0; i < gridSizeX; i++){
            for (int j = 0; j < gridSizeY; j++){
                neighbours[i, j] = new List<XpbdFluidParticel>();
            }
        }


        Vector3 s = new Vector3(-sizeX / 2, -sizeY / 2, 0);
        foreach (XpbdFluidParticel particle in particles)
        {
            int gridX = Mathf.FloorToInt((particle.x.x - s.x) / smoothingRadius);
            int gridY = Mathf.FloorToInt((particle.x.y - s.y) / smoothingRadius);

            gridX = Mathf.Clamp(gridX, 0, gridSizeX - 1);
            gridY = Mathf.Clamp(gridY, 0, gridSizeY - 1);

            neighbours[gridX, gridY].Add(particle);
        }
    }
    private List<XpbdFluidParticel> getNeighbours(Vector3 position)
    {
        Vector3 s = new Vector3(-sizeX / 2, -sizeY / 2, 0);
        int gridX = Mathf.FloorToInt((position.x - s.x) / smoothingRadius);
        int gridY = Mathf.FloorToInt((position.y - s.y) / smoothingRadius);

        List<XpbdFluidParticel> neighborParticles = new List<XpbdFluidParticel>();

        for (int x = -1; x <= 1; x++)
        {
            for (int y = -1; y <= 1; y++)
            {
                int neighborX = gridX + x;
                int neighborY = gridY + y;

                if (neighborX >= 0 && neighborX < gridSizeX && neighborY >= 0 && neighborY < gridSizeY)
                {
                    neighborParticles.AddRange(neighbours[neighborX, neighborY]);
                }
            }
        }

        return neighborParticles;
    }


    //---------------------------------------------------------------------------------------------//
    //-------------------------------- Rendering --------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    private void renderParicles()
    {
        Matrix4x4[] positions = new Matrix4x4[particles.Count];
        for (int i = 0; i < particles.Count; i++){
            positions[i] = Matrix4x4.TRS(particles[i].x, Quaternion.identity, Vector3.one);
        }

        Graphics.DrawMeshInstanced(particleMesh, 0, particleMaterial, positions);

    }
    private void OnDrawGizmos()
    {
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(Vector3.zero, new Vector3(sizeX, sizeY, 0));
        //Gizmos.DrawWireSphere(mousePos, mouseRadius);
    }

    //---------------------------------------------------------------------------------------------//
    //-------------------------------- Interaction ------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    Vector3 InteractionForce(Vector3 inputPos, float radius, float strength, int particleIndex)
    {
        Vector3 interactionForce = Vector2.zero;
        Vector3 offset = inputPos - particles[particleIndex].x;
        float sqrDst = Vector2.Dot(offset, offset);

        if (sqrDst < radius * radius)
        {
            float dst = Mathf.Sqrt(sqrDst);
            Vector3 dirToInputPoint = dst <= float.Epsilon ? Vector2.zero : offset / dst;
            float centreT = 1 - dst / radius;
            interactionForce += (dirToInputPoint * strength - particles[particleIndex].v) * centreT;
        }

        return interactionForce;
    }
    void mouseInput()
    {
        Vector3 mouseScreenPosition = Input.mousePosition;
        mouseScreenPosition.z = 20;
        mousePos = Camera.main.ScreenToWorldPoint(mouseScreenPosition);

        if (Input.GetMouseButtonDown(0)) strength = mouseStrength;

        else if (Input.GetMouseButtonDown(1)) strength = -mouseStrength;

        else if (Input.GetMouseButtonUp(0) || Input.GetMouseButtonUp(1)) strength = 0;
    }

    

}

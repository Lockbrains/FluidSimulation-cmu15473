using System.Runtime.InteropServices;
using UnityEngine;
using System.Collections.Generic;
using UnityEngine.Rendering;

[System.Serializable]


public class SPH3D : MonoBehaviour
{
    public struct Particle 
    {
        public float pressure; // 4
        public float density; // 8
        public Vector3 currentForce; // 20
        public Vector3 velocity; // 32
        public Vector3 position; // 44
    }

    [Header("SPH Constants")]
    public float deltaT = 0.01f;
    public float influenceRadius = 0.1f;
    private float pi = Mathf.PI;
    private float G = 9.81f;    

    [Header("Spawn Data")]
    public static Vector3Int numToSpawn = new Vector3Int(2,2,2);
    public Vector3 boundingBox = new Vector3(30,30,30);
    public Vector3 spawnBoxCenter = new Vector3(0,3,0);
    public Vector3 spawnBox = new Vector3(10,10,10);
    public float particleRadius = 1f;

    [Header("Rendering Data")]
    public Mesh mesh;
	public Material material;

    [Header("Compute Shader Data")]
    public ComputeShader computeShader;

    [Header("Particle Data")]

    public Vector3[] particles;
    public Vector3[] velocities;
    public Vector3[] forces;
    public float[] pressures;
    public float[] densities;

    private int particlesCount = numToSpawn.x * numToSpawn.y * numToSpawn.z;
    private int threadsPerGroup = 8;


    private ComputeBuffer particlesBuffer;
    private ComputeBuffer velocitiesBuffer;
    private ComputeBuffer forcesBuffer;
    private ComputeBuffer pressuresBuffer;
    private ComputeBuffer densitiesBuffer;
    private ComputeBuffer argsBuffer;

    

    private int integrateID;
    private int computeDensityID;
    private int computeForcesID; 
    




    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Awake()
    {
        SpawnParticlesInBox();
        BuildComputeBuffers();

        material.SetFloat("ParticleRadius", particleRadius);
        material.SetBuffer("Particles", particlesBuffer);

        
    }

    
    private void BuildComputeBuffers()
    {
        uint[] args = { mesh.GetIndexCount(0), (uint)particlesCount, mesh.GetIndexStart(0), mesh.GetBaseVertex(0), 0 };
        argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        argsBuffer.SetData(args);

        particlesBuffer = new ComputeBuffer(particlesCount, sizeof(float) * 3);
        particlesBuffer.SetData(particles);

        velocitiesBuffer = new ComputeBuffer(particlesCount, sizeof(float) * 3);
        velocitiesBuffer.SetData(velocities);

        forcesBuffer = new ComputeBuffer(particlesCount, sizeof(float) * 3);
        forcesBuffer.SetData(forces);

        pressuresBuffer = new ComputeBuffer(particlesCount, sizeof(float));
        pressuresBuffer.SetData(pressures);

        densitiesBuffer = new ComputeBuffer(particlesCount, sizeof(float));
        densitiesBuffer.SetData(densities);

        computeShader.SetFloat("DeltaTime", deltaT);
        computeShader.SetFloat("InfluenceRadius", influenceRadius);
        computeShader.SetFloat("Pi", pi);
        computeShader.SetFloat("G", G);

        integrateID = computeShader.FindKernel("Integrate");
        computeForcesID = computeShader.FindKernel("ComputeForces");
        computeDensityID = computeShader.FindKernel("ComputeDensityAndPressure");
        
        computeShader.SetBuffer(integrateID, "particles", particlesBuffer);
        computeShader.SetBuffer(integrateID, "velocities", velocitiesBuffer);
        computeShader.SetBuffer(integrateID, "forces", forcesBuffer);
        computeShader.SetBuffer(integrateID, "pressures", pressuresBuffer);
        computeShader.SetBuffer(integrateID, "densities", densitiesBuffer);

        computeShader.SetBuffer(computeForcesID, "particles", particlesBuffer);
        computeShader.SetBuffer(computeForcesID, "velocities", velocitiesBuffer);
        computeShader.SetBuffer(computeForcesID, "forces", forcesBuffer);
        computeShader.SetBuffer(computeForcesID, "pressures", pressuresBuffer);
        computeShader.SetBuffer(computeForcesID, "densities", densitiesBuffer);

        computeShader.SetBuffer(computeDensityID, "particles", particlesBuffer);
        computeShader.SetBuffer(computeDensityID, "velocities", velocitiesBuffer);
        computeShader.SetBuffer(computeDensityID, "forces", forcesBuffer);
        computeShader.SetBuffer(computeDensityID, "pressures", pressuresBuffer);
        computeShader.SetBuffer(computeDensityID, "densities", densitiesBuffer);


        

    }

    // Update is called once per frame
    private void Update(){
        material.SetBuffer("Particles", particlesBuffer);
		Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(500.0f, 500.0f, 500.0f)), argsBuffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
	}
    
    private void FixedUpdate()
    {   
        
        computeShader.Dispatch(computeForcesID, particlesCount/threadsPerGroup, 1, 1);
        computeShader.Dispatch(integrateID, particlesCount/threadsPerGroup, 1, 1);
        particlesBuffer.GetData(particles);

        /**
        computeForces();
        integrate();
        particlesBuffer.SetData(particles);
        **/

        
    }
    
    
    private void SpawnParticlesInBox()
    {
        Vector3 spawnTopLeft = spawnBoxCenter - spawnBox / 2;
        particles = new Vector3[particlesCount];
        velocities = new Vector3[particlesCount];
        forces = new Vector3[particlesCount];
        pressures = new float[particlesCount];
        densities = new float[particlesCount];

        for (int x = 0; x < numToSpawn.x; x++)
        {
            for (int y = 0; y < numToSpawn.y; y++)
            {
                for (int z = 0; z < numToSpawn.z; z++)
                {
                    Vector3 spawnPosition = spawnTopLeft + new Vector3(x * particleRadius * 2, y * particleRadius * 2, z * particleRadius * 2);
                    int i = x + y * numToSpawn.x + z * numToSpawn.x * numToSpawn.y;
                    particles[i] = spawnPosition;
                    velocities[i] = Vector3.zero;
                    forces[i] = Vector3.zero;
                    pressures[i] = 0;
                    densities[i] = 0;
                
                }
            }
        }
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.blue;
        Gizmos.DrawWireCube(Vector3.zero, boundingBox);
        
        Gizmos.color = Color.cyan;
        Gizmos.DrawWireCube(spawnBoxCenter, spawnBox);
        

    }

    private void integrate(){
        for (int i = 0; i < particlesCount; i++)
        {
        
            velocities[i] += forces[i] * deltaT;
            particles[i] += velocities[i] * deltaT;
            forces[i] = Vector3.zero;
        }
    }
    private void computeForces(){
        for (int i = 0; i < particlesCount; i++)
        {
            forces[i] = G * Vector3.down;
        }
    }
    
}

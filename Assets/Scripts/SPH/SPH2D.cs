using System.Runtime.InteropServices;
using UnityEngine;
using System.Collections.Generic;
using UnityEngine.Rendering;
using System;

[System.Serializable]
[StructLayout(LayoutKind.Sequential, Size = 36)]
public struct Particle
{
    public float pressure; // 4
    public float density; // 8
    public Vector2 currentForce; // 16
    public Vector2 velocity; // 24
    public Vector2 position; // 32
    public int cell; // 36
}

public class SPH2D : MonoBehaviour
{

    [Header("SPH Constants")]
    public float mass = 1.0f;
    public float deltaT = 0.005f;
    public float bounce = 0.8f;
    private float pi = Mathf.PI;
    private float G = 9.81f;
    public float k = 1;
    public float restDensity = 0.1f;   
    public float viscosity = 0.05f;
     

    public float influenceRadius = 3.0f;

    [Header("Spawn Data")]
    public Vector2Int numToSpawn = new Vector2Int(128,96);
    public Vector2 boundingBox = new Vector2(300,200);
    public Vector2 spawnBoxCenter = new Vector2(0,0);
    public Vector2 spawnBox = new Vector2(256,192);
    public float particleRadius = 1f;

    [Header("Rendering Data")]
    public Mesh mesh;
	public Material material;

    [Header("Compute Shader Data")]
    public ComputeShader computeShader;

    [Header("Particle Data")]

    public Particle[] particles;
    public int[] cellStart;
    public int numCells = 10;

    private int particlesCount = 0;
    private int threadsPerGroup = 128;


    public ComputeBuffer particlesBuffer;
    public ComputeBuffer cellStartBuffer;
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

        particlesBuffer = new ComputeBuffer(particlesCount, Marshal.SizeOf(typeof(Particle)));
        particlesBuffer.SetData(particles);

        cellStartBuffer = new ComputeBuffer(numCells * numCells, Marshal.SizeOf(typeof(int)));
        cellStartBuffer.SetData(cellStart);


        computeShader.SetFloat("Mass", mass);
        computeShader.SetFloat("DeltaTime", deltaT);
        computeShader.SetFloat("Bounce", bounce);
        computeShader.SetFloat("Pi", pi);
        computeShader.SetFloat("G", G);
        computeShader.SetFloat("K", k);
        computeShader.SetFloat("RestDensity", restDensity);
        computeShader.SetFloat("Viscosity", viscosity);
        computeShader.SetInt("NumCells", numCells);


        Vector4 boundingBox = new Vector4(this.boundingBox.x, this.boundingBox.y, 0, 0);
        computeShader.SetVector("BoundingBox", boundingBox);


        float influenceRadius2 = influenceRadius * influenceRadius;
        float influenceRadius3 = influenceRadius * influenceRadius2;
        float influenceRadius6 = influenceRadius3 * influenceRadius3;
        float influenceRadius9 = influenceRadius6 * influenceRadius3;
        computeShader.SetFloat("InflRad", influenceRadius);
        computeShader.SetFloat("InflRad2", influenceRadius2);
        computeShader.SetFloat("InflRad3", influenceRadius3);
        computeShader.SetFloat("InflRad6", influenceRadius6);
        computeShader.SetFloat("InflRad9", influenceRadius9);

        integrateID = computeShader.FindKernel("Integrate");
        computeForcesID = computeShader.FindKernel("ComputeForces");
        computeDensityID = computeShader.FindKernel("ComputeDensityAndPressure");
        
        computeShader.SetBuffer(integrateID, "particles", particlesBuffer);

        computeShader.SetBuffer(computeForcesID, "particles", particlesBuffer);
        computeShader.SetBuffer(computeForcesID, "cellStart", cellStartBuffer);

        computeShader.SetBuffer(computeDensityID, "particles", particlesBuffer);
        computeShader.SetBuffer(computeDensityID, "cellStart", cellStartBuffer);



        

    }

    // Update is called once per frame
    private void Update(){
		Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(500.0f, 500.0f, 500.0f)), argsBuffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
         computeShader.Dispatch(computeDensityID, particlesCount/threadsPerGroup, 1, 1);
        computeShader.Dispatch(computeForcesID, particlesCount/threadsPerGroup, 1, 1);
        computeShader.Dispatch(integrateID, particlesCount/threadsPerGroup, 1, 1);

        computeShader.SetFloat("Mass", mass);
        computeShader.SetFloat("DeltaTime", deltaT);
        computeShader.SetFloat("Bounce", bounce);
        computeShader.SetFloat("Pi", pi);
        computeShader.SetFloat("G", G);
        computeShader.SetFloat("K", k);
        computeShader.SetFloat("RestDensity", restDensity);
        computeShader.SetFloat("Viscosity", viscosity);

        float influenceRadius2 = influenceRadius * influenceRadius;
        float influenceRadius3 = influenceRadius * influenceRadius2;
        float influenceRadius6 = influenceRadius3 * influenceRadius3;
        float influenceRadius9 = influenceRadius6 * influenceRadius3;
        computeShader.SetFloat("InflRad", influenceRadius);
        computeShader.SetFloat("InflRad2", influenceRadius2);
        computeShader.SetFloat("InflRad3", influenceRadius3);
        computeShader.SetFloat("InflRad6", influenceRadius6);
        computeShader.SetFloat("InflRad9", influenceRadius9);
        particlesBuffer.GetData(particles);
        Array.Sort(particles, (p1, p2) => p1.cell.CompareTo(p2.cell));
        for(int i = 0; i < particles.Length; i++){
            if(i == 0 || particles[i].cell != particles[i-1].cell){
                cellStart[particles[i].cell] = i;
            }
        }
        particlesBuffer.SetData(particles);
        cellStartBuffer.SetData(cellStart);

    }
    
    private void FixedUpdate()
    {   
       

    

        /**
        computeForces();
        integrate();
        particlesBuffer.SetData(particles);
        **/

        
    }
    
    
    private void SpawnParticlesInBox()
    {
        Vector2 spawnTopLeft = spawnBoxCenter - spawnBox / 2;
        List<Particle> particlesList = new List<Particle>();
        cellStart = new int[numCells * numCells];
        
        for (int x = 0; x < numToSpawn.x; x++)
        {
            for (int y = 0; y < numToSpawn.y; y++)
            {
                Vector2 spawnPosition = spawnTopLeft + new Vector2(x * particleRadius * 2, y * particleRadius * 2);
                spawnPosition += new Vector2(UnityEngine.Random.Range(-particleRadius, particleRadius), 
                UnityEngine.Random.Range(-particleRadius, particleRadius));
                int posCell = (int)(((spawnPosition.x+boundingBox.x/2)/(boundingBox.x/numCells)) + ((int)(spawnPosition.y+boundingBox.y/2)/(boundingBox.y/numCells))*numCells);
                Particle p  = new Particle{position = spawnPosition, cell = posCell}; 
                particlesList.Add(p); 
            }
        }
        particlesList.Sort((p1, p2) => p1.cell.CompareTo(p2.cell));
        for(int i = 0; i < particlesList.Count; i++){
            if(i == 0 || particlesList[i].cell != particlesList[i-1].cell){
                cellStart[particlesList[i].cell] = i;
            }
        }
        particles = particlesList.ToArray();
        particlesCount = particles.Length;
    }
    
    private void OnDrawGizmos()
    {
        Gizmos.color = Color.blue;
        Vector3 box = new Vector3(boundingBox.x, boundingBox.y, 0);
        Gizmos.DrawWireCube(Vector3.zero, box);
        
        if(!Application.isPlaying){
            Gizmos.color = Color.cyan;
            Vector3 sBox = new Vector3(spawnBox.x, spawnBox.y, 0);
            Vector3 sBoxCenter = new Vector3(spawnBoxCenter.x, spawnBoxCenter.y, 0);
            Gizmos.DrawWireCube(sBoxCenter, sBox);
        }
        
    }
/**
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
    */
}

using UnityEngine;

public class Eulerian2D : MonoBehaviour
{
    public float[] verticalVelocities;
    private float[] copyVert;
    public float[] horizontalVelocities;
    private float[] copyHor;
    public int[] vertBound;
    public int[] horBound;
    public float[] densities;
    private float[] copyDens;
    public Vector2[] pos;

    [Header("Params")]
    public int width = 10;
    public int height = 10;
    public float deltaT = 0.005f;
    public float maxDensity = 100.0f;
    public float cellSize = 1.0f;
    public float overRelaxation = 1.9f;
    public int divIter = 100;

    public ComputeShader computeShader;

    private ComputeBuffer vertVelBuffer;
    private ComputeBuffer horVelBuffer;
    private ComputeBuffer vertBoundBuffer;
    private ComputeBuffer horBoundBuffer;
    private ComputeBuffer densityBuffer;
    private ComputeBuffer positionBuffer;
    public ComputeBuffer argsBuffer;

    [Header("Rendering Data")]
    public Mesh mesh;
	public Material material;

    private  bool[] clicked;

    private int numThreads = 10;
    private int kernelAdvect;
    private int kernelAdvectDen;
    private int kernelClearDiv1;
    private int kernelClearDiv2;



    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Awake()
    {
        verticalVelocities = new float[width * (height+1)];
        copyVert = new float[width * (height+1)];
        horizontalVelocities = new float[(width+1) * height];
        copyHor = new float[(width+1) * height];
        vertBound = new int[width * (height+1)];
        horBound = new int[(width+1) * height];
        densities = new float[width * height];
        copyDens = new float[width * height];
        pos = new Vector2[width * height];
        clicked = new bool[width * height];
        for(int i=0; i<height+1; i++){
            for(int j=0; j<width; j++){
                verticalVelocities[j + i * width] = 0;
                copyVert[j + i * width] = 0;
                if(i==0 || i==height){
                    vertBound[j + i * width] = 1;
                }
                else{
                    vertBound[j + i * width] = 0;
                }
            }
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width+1; j++){
                horizontalVelocities[j + i * (width+1)] = 0;
                copyHor[j + i * (width+1)] = 0;
                if(j==0 || j==width){
                    horBound[j + i * (width+1)] = 1;
                }
                else{
                    horBound[j + i * (width+1)] = 0;
                }

            }
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                densities[j + i * width] = 0;
                copyDens[j + i * width] = 0;
                pos[j + i * width] = new Vector2(j, i);
            }
        }

        buildBuffers();
    }

void buildBuffers(){
     // Create buffers
        vertVelBuffer = new ComputeBuffer(width * (height+1), sizeof(float));
        vertVelBuffer.SetData(verticalVelocities);
        horVelBuffer = new ComputeBuffer((width+1) * height, sizeof(float));
        horVelBuffer.SetData(horizontalVelocities);
        vertBoundBuffer = new ComputeBuffer(width * (height+1), sizeof(int));
        vertBoundBuffer.SetData(vertBound);
        horBoundBuffer = new ComputeBuffer((width+1) * height, sizeof(int));
        horBoundBuffer.SetData(horBound);
        densityBuffer = new ComputeBuffer(width * height, sizeof(float));
        densityBuffer.SetData(densities);
        positionBuffer = new ComputeBuffer(width * height, sizeof(float) * 2);
        positionBuffer.SetData(pos);

        // Set buffers in material
        uint[] args = { mesh.GetIndexCount(0), (uint)(width * height), mesh.GetIndexStart(0), mesh.GetBaseVertex(0), 0};
        argsBuffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        argsBuffer.SetData(args);
        material.SetBuffer("DensityBuffer", densityBuffer);
        material.SetBuffer("PositionBuffer", positionBuffer);
        material.SetFloat("_MaxDensity", maxDensity);
        material.SetFloat("_CellSize", cellSize);

        // Find kernel IDs
        kernelAdvect = computeShader.FindKernel("Advect");
        kernelAdvectDen = computeShader.FindKernel("AdvectDensity");
        kernelClearDiv1 = computeShader.FindKernel("ClearDivergence1");
        kernelClearDiv2 = computeShader.FindKernel("ClearDivergence2");

        // Set buffers in compute shader
        computeShader.SetBuffer(kernelClearDiv1, "vertVel", vertVelBuffer);
        computeShader.SetBuffer(kernelClearDiv1, "horizVel", horVelBuffer);
        computeShader.SetBuffer(kernelClearDiv1, "vertBoundary", vertBoundBuffer);
        computeShader.SetBuffer(kernelClearDiv1, "horizBoundary", horBoundBuffer);

        computeShader.SetBuffer(kernelClearDiv2, "vertVel", vertVelBuffer);
        computeShader.SetBuffer(kernelClearDiv2, "horizVel", horVelBuffer);
        computeShader.SetBuffer(kernelClearDiv2, "vertBoundary", vertBoundBuffer);
        computeShader.SetBuffer(kernelClearDiv2, "horizBoundary", horBoundBuffer);

        // Set variables in compute shader
        computeShader.SetInt("Width", width);
        computeShader.SetInt("Height", height);
        computeShader.SetFloat("DeltaT", deltaT);
        computeShader.SetFloat("MaxDensity", maxDensity);
        computeShader.SetFloat("CellSize", cellSize);
        computeShader.SetFloat("OverRelaxation", overRelaxation);

}

void opaqueCell(int x, int y)
{
    densities[x + y * width] = 0;
    vertBound[x + y * width] = 1;
    vertBound[x + (y+1) * width] = 1;
    horBound[x + y * (width+1)] = 1;
    horBound[x+1 + y * (width+1)] = 1;
    verticalVelocities[x + y * width] = 0;
    verticalVelocities[x + (y+1) * width] = 0;
    horizontalVelocities[x + y * (width+1)] = 0;
    horizontalVelocities[x+1 + y * (width+1)] = 0;

}
void HandleMouseClick(int button)
{
    Vector3 mouseScreenPos = Input.mousePosition;
    Vector3 mouseWorldPos = Camera.main.ScreenToWorldPoint(new Vector3(mouseScreenPos.x, mouseScreenPos.y, Camera.main.nearClipPlane));
    mouseWorldPos.z = 0; // Ensure it's in 2D space

    // Convert world position to grid indices
    int cellX = Mathf.FloorToInt(mouseWorldPos.x);
    int cellY = Mathf.FloorToInt(mouseWorldPos.y);

    // Check if within bounds
    if (cellX >= 0 && cellX < width && cellY >= 0 && cellY < height)
    {
        if (button == 0)
        {
            //IncreaseVelocityAtCell(cellX, cellY, false);
            opaqueCell(cellX, cellY);
            clicked[cellX + cellY * width] = !clicked[cellX + cellY * width];

        }
        else if (button == 1)
        {
            IncreaseVelocityAtCell(cellX, cellY, true);
            clicked[cellX + cellY * width] = !clicked[cellX + cellY * width];
        }
    }
}

    // Update is called once per frame


    void Update()
    {
        if (Input.GetMouseButtonDown(0))
        {
            HandleMouseClick(0);
        }
        else if (Input.GetMouseButtonDown(1))
        {
            HandleMouseClick(1);
        }

        Graphics.DrawMeshInstancedIndirect(mesh, 0, material, new Bounds(Vector3.zero, new Vector3(width, height, 0)), argsBuffer, castShadows: UnityEngine.Rendering.ShadowCastingMode.Off);
        
    }


    void FixedUpdate()
    {
        
        horBound[height/2 * (width+1)] = 0;
        horBound[height/2 * (width+1) - 1] = 0;
        densities[height/2 * width] = maxDensity;
        horizontalVelocities[height/2 * (width+1)] = 5.0f;
        advectDensity();
        advect();
        vertVelBuffer.SetData(verticalVelocities);
        horVelBuffer.SetData(horizontalVelocities);
        
        for(int n=0; n<divIter; n++){
            computeShader.Dispatch(kernelClearDiv1, width*height/numThreads, 1, 1);
            computeShader.Dispatch(kernelClearDiv2, width*height/numThreads, 1, 1);
        }
        
        
        densityBuffer.SetData(densities);
        material.SetBuffer("DensityBuffer", densityBuffer);

        vertVelBuffer.GetData(verticalVelocities);
        horVelBuffer.GetData(horizontalVelocities);

            
           
        
    }

    
    void clearDivergence()
    {   
        for(int n=0; n<divIter; n++){
            for(int i=0; i<height; i++){
                for(int j=0; j<width; j++){
                    float div = -verticalVelocities[j+i*width]+verticalVelocities[j+(i+1)*width] + 
                                -horizontalVelocities[j+i*(width+1)]+horizontalVelocities[j+1+i*(width+1)];
                    div *= overRelaxation;
                    int vup = vertBound[j + (i+1) * width]==1 ? 0 : 1;
                    int vdown = vertBound[j + i * width]==1 ? 0 : 1;
                    int hleft = horBound[j + i * (width+1)]==1 ? 0 : 1;
                    int hright = horBound[j+1 + i * (width+1)]==1 ? 0 : 1;
                    int tot = vup+vdown+hleft+hright;
                    if(tot == 0) continue;
                    


                    verticalVelocities[j+i*width]+=vdown*div/tot;
                    verticalVelocities[j+(i+1)*width]-=vup*div/tot;
                    horizontalVelocities[j+i*(width+1)]+=hleft*div/tot;
                    horizontalVelocities[j+1+i*(width+1)]-=hright*div/tot;
                }
            }
        }
    }

    float avgVert(int x, int y)
    {
        // up left, up right etc.
        int yul = (x-1) >= 0 ? (x-1) + (y+1)*width : -1;
        int yur = x < width ? x + (y+1)*width : -1;
        int ydl = (x-1) >= 0 ? (x-1) + y*width : -1;
        int ydr = x < width ? x + y*width : -1;

        float vel = 0;
        int count = 0;
        if(yul != -1){
            vel += verticalVelocities[yul];
            count++;
        }
        if(yur != -1){
            vel += verticalVelocities[yur];
            count++;
        }
        if(ydl != -1){
            vel += verticalVelocities[ydl];
            count++;
        }
        if(ydr != -1){
            vel += verticalVelocities[ydr];
            count++;
        }
        return vel/count;
        
    }

     float avgHor(int x, int y)
    {
        // up left, up right etc.
        int xul = y < height ? x + y*(width+1) : -1;
        int xur = y < height ? x+1 + y*(width+1) : -1;
        int xdl = (y-1) >= 0 ? x + (y-1)*(width+1) : -1;
        int xdr = (y-1) >= 0 ? x+1 + (y-1)*(width+1) : -1;

        float vel = 0;
        int count = 0;
        if(xul != -1){
            vel += horizontalVelocities[xul];
            count++;
        }
        if(xur != -1){
            vel += horizontalVelocities[xur];
            count++;
        }
        if(xdl != -1){
            vel += horizontalVelocities[xdl];
            count++;
        }
        if(xdr != -1){
            vel += horizontalVelocities[xdr];
            count++;
        }
        return vel/count;
    }
    float bilineraInter(float x, float y, float[] data, int width, int height)
    {
        int x0 = Mathf.FloorToInt(x);
        int x1 = Mathf.Min(x0 + 1, width - 1);
        int y0 = Mathf.FloorToInt(y);
        int y1 = Mathf.Min(y0 + 1, height - 1);

        float wx = x-x0;
        float wy = y-y0;

        // Bilinear interpolation
        return 
            (1 - wx) * (1 - wy) * data[x0 + y0 * width] +
            wx * (1 - wy) * data[x1 + y0 * width] +
            (1 - wx) * wy * data[x0 + y1 * width] +
            wx * wy * data[x1 + y1 * width];
    }
    void advect()
    {   
        //vertical
        for(int i=0; i<height+1; i++){
            for(int j=0; j<width; j++){
                // In the middle of the edge x axis
                if(vertBound[j + i * width] == 1){
                    copyVert[j + i * width] = 0;
                    continue;
                }
                float x = j + 0.5f - deltaT*avgHor(j, i);
                float y = i - deltaT*verticalVelocities[j+i*width];

                x = Mathf.Clamp(x, 0, width-1);
                y = Mathf.Clamp(y, 0, height);

               

                // Bilinear interpolation
                copyVert[j + i * width] = bilineraInter(x, y, verticalVelocities, width, height+1);
                   
            }
        }

        //horizontal
        for(int i=0; i<height; i++){
            for(int j=0; j<width+1; j++){
                // In the middle of the edge y axis
                if(horBound[j + i * (width+1)] == 1){
                    copyHor[j + i * (width+1)] = 0;
                    continue;
                }
                float x = j - deltaT*horizontalVelocities[j+i*(width+1)];
                float y = i+ 0.5f - deltaT*avgVert(j, i);

                x = Mathf.Clamp(x, 0, width);
                y = Mathf.Clamp(y, 0, height-1);

               

                // Bilinear interpolation
                copyHor[j + i * (width + 1)] = bilineraInter(x, y, horizontalVelocities, width+1, height);
                
            }
        }
        // Swap buffers
        var temp = verticalVelocities;
        verticalVelocities = copyVert;
        copyVert = temp;

        temp = horizontalVelocities;
        horizontalVelocities = copyHor;
        copyHor = temp;
        
    }



    public float res;
    void advectDensity()
    {
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                float x = j - deltaT*(horizontalVelocities[j+i*(width+1)] + horizontalVelocities[j+1+i*(width+1)])/2;
                float y = i - deltaT*(verticalVelocities[j+i*width] + verticalVelocities[j+(i+1)*width])/2;

                x = Mathf.Clamp(x, 0, width-1);
                y = Mathf.Clamp(y, 0, height-1);

               

                // Bilinear interpolation
                res = bilineraInter(x, y, densities, width, height);
                copyDens[j + i * width] = res;
            }
        }

        // Swap buffers
        var temp = densities;
        densities = copyDens;
        copyDens = temp;
    }

    void IncreaseVelocityAtCell(int x, int y, bool vert)
    {
        if(vert){
            verticalVelocities[x + y * width] += 5.0f;
        }
        else{
            horizontalVelocities[x + y * (width + 1)] += 5.0f;
        }
       
    }

    void OnDrawGizmos()
    {
        if (verticalVelocities == null || horizontalVelocities == null) return;
        if(Application.isPlaying){
            
            Gizmos.color = Color.white;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {   
                    /*
                    float div = verticalVelocities[j+i*width]-verticalVelocities[j+(i+1)*width] + 
                                horizontalVelocities[j+i*(width+1)]-horizontalVelocities[j+1+i*(width+1)];
                    if(div < -0.001 || div > 0.001){
                        Gizmos.color = Color.red;
                        Gizmos.DrawWireCube(new Vector3(j , i , 0), new Vector3(1, 1, 0));
                        Gizmos.color = Color.white;
                    }
                    */
                    // Compute position at the center of the cell
                    Vector3 position = new Vector3(j + 0.5f, i + 0.5f, 0);

                    // Compute the average velocity at the cell center
                    float averageVertical = (verticalVelocities[j + i * width]
                                        + verticalVelocities[j + (i + 1) * width]) * 0.5f;

                    float averageHorizontal = (horizontalVelocities[j + i * (width + 1)]
                                            + horizontalVelocities[j + 1 + i * (width + 1)]) * 0.5f;

                    Vector3 velocity = new Vector3(averageHorizontal, averageVertical, 0);

                    // Draw the velocity vector as an arrow
                    Vector3 arrowEnd = position + velocity * 0.5f; // Scale velocity for visibility
                    Gizmos.DrawLine(position, arrowEnd);

                    // Draw arrowhead
                    Vector3 arrowDirection = velocity.normalized * 0.1f; // Normalize for arrowhead size
                    Vector3 perpendicular = new Vector3(-arrowDirection.y, arrowDirection.x, 0); // Perpendicular vector

                    Gizmos.DrawLine(arrowEnd, arrowEnd - arrowDirection + perpendicular * 0.5f);
                    Gizmos.DrawLine(arrowEnd, arrowEnd - arrowDirection - perpendicular * 0.5f);

                    // Optionally draw a square for each cell
                    Gizmos.color = Color.gray;
                    Gizmos.DrawWireCube(new Vector3(j + 0.5f, i + 0.5f, 0), new Vector3(1, 1, 0));
                        
                    Gizmos.color = Color.white; // Reset color for arrows
                }
                    
            }
        }
        else{
            for(int i=0; i<height; i++){
                for(int j=0; j<width; j++){
                    // Optionally draw a square for each cell
                    Gizmos.color = Color.gray;
                    Gizmos.DrawWireCube(new Vector3(j + 0.5f, i + 0.5f, 0), new Vector3(1, 1, 0));
                    Gizmos.color = Color.blue; // Reset color for arrows
                    
                }
            }
        }
    }

    
}

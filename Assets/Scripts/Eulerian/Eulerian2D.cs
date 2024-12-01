using UnityEngine;

public class Eulerian2D : MonoBehaviour
{
    public float[] verticalVelocities;
    public float[] copyVert;
    public float[] horizontalVelocities;
    public float[] copyHor;

    public uint width = 10;
    public uint height = 10;
    public float deltaT = 0.005f;
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        verticalVelocities = new float[width * (height+1)];
        copyVert = new float[width * (height+1)];
        horizontalVelocities = new float[(width+1) * height];
        copyHor = new float[(width+1) * height];

        // Example: Initialize with random small velocities
        for (int i = 0; i < width * height; i++)
        {
            verticalVelocities[i] = Random.Range(-0.1f, 0.1f);
            horizontalVelocities[i] = Random.Range(-0.1f, 0.1f); 
        }
    }

    // Update is called once per frame
    void Update()
    {
        
        
    }

    void FixedUpdate()
    {
        advect();
        clearDivergence();
    }

    void clearDivergence()
    {
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                float div = verticalVelocities[j+i*width]-verticalVelocities[j+(i+1)*width] + 
                            horizontalVelocities[j+i*width]-horizontalVelocities[j+1+i*width];
                verticalVelocities[j+i*width]-=div/4;
                verticalVelocities[j+(i+1)*width]+=div/4;
                horizontalVelocities[j+i*width]-=div/4;
                horizontalVelocities[j+1+i*width]+=div/4;
            }
        }
    }

    float avgVert(int x, int y)
    {
        // up left, up right etc.
        int yul = (int)(Mathf.Clamp(x-1, 0, width-1) + y*width);
        int yur = (int)(Mathf.Clamp(x, 0, width-1) + y*width);
        int ydl = (int)(Mathf.Clamp(x-1, 0, width-1) + (y+1)*width);
        int ydr = (int)(Mathf.Clamp(x, 0, width-1) + (y+1)*width);
        return (verticalVelocities[yul]+verticalVelocities[yur]+verticalVelocities[ydl]+verticalVelocities[ydr])/4;
    }

     float avgHor(int x, int y)
    {
        // up left, up right etc.
        int xul = (int)(x + Mathf.Clamp(y-1, 0, height-1)*(width+1));
        int xur = (int)(x+1 + Mathf.Clamp(y-1, 0, height-1)*(width+1));
        int xdl = (int)(x + Mathf.Clamp(y, 0, height-1)*(width+1));
        int xdr = (int)(x+1 + Mathf.Clamp(y, 0, height-1)*(width+1));
        return (horizontalVelocities[xul]+horizontalVelocities[xur]+horizontalVelocities[xdl]+horizontalVelocities[xdr])/4;
    }
    void advect()
    {   
        //vertical
        for(int i=0; i<height+1; i++){
            for(int j=0; j<width; j++){
                // In the middle of the edge x axis
                float x = j + 0.5f - deltaT*avgHor(j, i);
                float y = i - deltaT*verticalVelocities[j+i*width];

                x = Mathf.Clamp(x, 0, width-1);
                y = Mathf.Clamp(y, 0, height);

                int x0 = Mathf.FloorToInt(x);
                int x1 = x0+1;
                x1 = (int)Mathf.Clamp(x1, 0, width - 1);
                int y0 = Mathf.FloorToInt(y);
                int y1 = y0+1;
                y1 = (int)Mathf.Clamp(y1, 0, height);

                float wx = x-x0;
                float wy = y-y0;

                // Bilinear interpolation
                copyVert[j + i * width] = 
                    (1 - wx) * (1 - wy) * verticalVelocities[x0 + y0 * width] +
                    wx * (1 - wy) * verticalVelocities[x1 + y0 * width] +
                    (1 - wx) * wy * verticalVelocities[x0 + y1 * width] +
                    wx * wy * verticalVelocities[x1 + y1 * width];
            }
        }

        //horizontal
        for(int i=0; i<height; i++){
            for(int j=0; j<width+1; j++){
                // In the middle of the edge y axis
                float x = j - deltaT*horizontalVelocities[j+i*(width+1)];
                float y = i+ 0.5f - deltaT*avgVert(j, i);

                x = Mathf.Clamp(x, 0, width);
                y = Mathf.Clamp(y, 0, height-1);

                int x0 = Mathf.FloorToInt(x);
                int x1 = x0+1;
                x1 = (int)Mathf.Clamp(x1, 0, width);
                int y0 = Mathf.FloorToInt(y);
                int y1 = y0+1;
                y1 = (int)Mathf.Clamp(y1, 0, height-1);

                float wx = x-x0;
                float wy = y-y0;

                // Bilinear interpolation
                copyHor[j + i * (width + 1)] = 
                (1 - wx) * (1 - wy) * horizontalVelocities[x0 + y0 * (width + 1)] +
                wx * (1 - wy) * horizontalVelocities[x1 + y0 * (width + 1)] +
                (1 - wx) * wy * horizontalVelocities[x0 + y1 * (width + 1)] +
                wx * wy * horizontalVelocities[x1 + y1 * (width + 1)];
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

    void OnDrawGizmos()
    {
        if (verticalVelocities == null || horizontalVelocities == null) return;

        Gizmos.color = Color.blue;
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                // Compute position at the center of the cell
                Vector3 position = new Vector3(j + 0.5f, i + 0.5f, 0); // Center of the cell

                // Compute the average velocity at the cell center
                float averageVertical = (verticalVelocities[j + i * width] 
                                    + verticalVelocities[j + (i + 1) * width]) * 0.5f;

                float averageHorizontal = (horizontalVelocities[j + i * (width + 1)] 
                                        + horizontalVelocities[j + 1 + i * (width + 1)]) * 0.5f;

                Vector3 velocity = new Vector3(averageHorizontal, averageVertical, 0);

                // Draw the velocity vector
                Gizmos.DrawLine(position, position + velocity*100); // Scale for visibility
            }
        }
    }

    
}

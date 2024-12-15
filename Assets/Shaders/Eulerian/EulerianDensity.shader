// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Custom/Eulerian"
{
    Properties
    {
        _MaxDensity("Max Density", Float) = 100.0
        _CellSize("Cell Size", Float) = 1.0
    }

    SubShader
    {
        Tags { "RenderType" = "Opaque" }
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag


            #include "UnityCG.cginc"

            // Structured buffer for the positions of the grid cells
            StructuredBuffer<float2> PositionBuffer;

            //  Structured buffer for the density of the grid cells
            StructuredBuffer<float> DensityBuffer;
            float _MaxDensity;
            float _CellSize;

            struct appdata_t
            {
                float4 vertex : POSITION;
                
            };

            struct v2f
            {
                float4 pos : SV_POSITION;
                float4 color : COLOR;
                
            };

            // Vertex shader
            v2f vert(appdata_t v, uint instanceID: SV_InstanceID)
            {
                v2f o;
                float4x4 modelMatrix = {1, 0, 0, PositionBuffer[instanceID].x+_CellSize/2,
                                         0, 1, 0, PositionBuffer[instanceID].y+_CellSize/2,
                                         0, 0, 1, 0,
                                         0, 0, 0, 1};
                float4 pos = mul(modelMatrix, v.vertex);
                o.pos = UnityObjectToClipPos(pos);


                float density = log(DensityBuffer[instanceID] + 1.0)/log(_MaxDensity+1.0);
                
                half r = density > 0.5 ? 2.0 * (density - 0.5) : 0.0; // Red starts increasing after density > 0.5
                half g = density <= 0.5 ? 2.0 * density : 2.0 * (1.0 - density); // Green peaks at density = 0.5
                half b = density <= 0.5 ? 1.0 - 2.0 * density : 0.0; // Blue fades out as density increases beyond 0.5

                //float density = clamp(DensityBuffer[instanceID]/_MaxDensity, 0.0, 1.0);
                //o.color = half4(density, 0, 1 - density, 1);
                r *= density;
                g *= density;
                b *= density;
                o.color = half4(r, g, b, 1);
                return o;
            }

            // Fragment shader
            half4 frag(v2f i, uint id : SV_InstanceID) : SV_Target{
         
                return i.color;
            }
            ENDCG
        }
    }

    Fallback "Diffuse"
}
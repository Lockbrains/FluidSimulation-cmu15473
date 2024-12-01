Shader "CustomRenderTexture/EulerianTexture"
{
    Properties
    {
        Velocity ("Velocity Texture (RGBA)", 2D) = "white" {}
        _Color ("Color", Color) = (1, 0, 0, 1) // Set color for the velocity representation
        gridWidth ("Grid Width", Int) = 128 
        gridHeight ("Grid Height", Int) = 128

    }
    SubShader
    {
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            // Texture to store velocity (RGBA) and Color property
            sampler2D Velocity;
            float4 _Color;
            uint gridWidth;
            uint gridHeight; 

            struct appdata_t
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            // Vertex Shader
            v2f vert(appdata_t v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                // Scale UVs to map the texture to the grid size
                o.uv = v.uv * float2(gridWidth, gridHeight);  // gridWidth and gridHeight are the dimensions of your grid/texture
                return o;
            }

            // Fragment Shader
            half4 frag(v2f i) : SV_Target
            {
                // Sample the velocity texture, assumes RGBA stores velocity components
                // (x, y, z = velocity components, w = not used)
                float4 velocity = tex2D(Velocity, i.uv);

                // Calculate velocity magnitude from RGB (x, y, z) components
                float speed = length(velocity.rgb); // This calculates the speed from the velocity components

                // Color visualization based on speed
                // You can customize the color further, here we scale it with the chosen color
                float3 color = _Color.rgb * speed; // Apply color scaling based on the magnitude

                // Return the final color as a RGBA value with a constant alpha
                return half4(color, 1.0);
            }
            ENDCG
        }
    }
}

Shader "FluidSimulation/VertexDisplacement"
{
    Properties
    {
        _HeightMap("Height Map", 2D) = "white" {} 
        _DisplacementScale("Displacement Scale", Float) = 1.0
        _MainTex("Main Texture", 2D) = "white" {}
        _LightColor("Light Color", Color) = (1,1,1,1)
        _AmbientColor("Ambient Color", Color) = (0.2, 0.2, 0.2, 1) 
        _SpecularPower("Specular Power", Float) = 16.0
        _ChoppyFactor("Choppy Factor", Float) = 0.2
    }

    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 200

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"
            
            sampler2D _HeightMap;
            sampler2D _MainTex;
            float _DisplacementScale;
            float4 _LightColor;
            float4 _AmbientColor;
            float _ChoppyFactor;
            float _SpecularPower;
            
            float4x4 _ObjectToWorld;
            float4x4 _WorldToObject;
            float4x4 _UnityMatrixVP;

            struct appdata
            {
                float4 vertex : POSITION; 
                float2 uv : TEXCOORD0;  
            };

            struct v2f
            {
                float4 vertex : SV_POSITION;
                float3 worldPos: TEXCOORD1;
                float2 uv : TEXCOORD0;       
            };


            v2f vert(appdata v)
            {
                v2f o;
                o.uv = v.uv;
                
                float height = tex2Dlod(_HeightMap, float4(v.uv, 0, 0)).g;
                float dx = tex2Dlod(_HeightMap, float4(v.uv, 0, 0)).r;
                float dz = tex2Dlod(_HeightMap, float4(v.uv, 0, 0)).b;
                
                float3 newPos = v.vertex.xyz;
                newPos.y += height * _DisplacementScale;
                newPos.x += dx * _ChoppyFactor;
                newPos.z += dz * _ChoppyFactor;

                float3 worldPos = UnityObjectToWorldDir(newPos);
                
                o.vertex = UnityObjectToClipPos(newPos);
                o.worldPos = worldPos;

                return o;
            }
            
            float4 frag(v2f i) : SV_Target
            {
                float3 ambient = _AmbientColor.rgb;
                
                float3 lightDir = normalize(float3(0, 1, -1));
                
                float heightL = tex2D(_HeightMap, i.uv + float2(-0.01, 0)).g;
                float heightR = tex2D(_HeightMap, i.uv + float2(0.01, 0)).g;
                float heightD = tex2D(_HeightMap, i.uv + float2(0, -0.01)).g;
                float heightU = tex2D(_HeightMap, i.uv + float2(0, 0.01)).g;

                float3 normal = normalize(float3(heightL - heightR, 2.0, heightD - heightU));
                
                float diff = max(0.0, dot(normal, lightDir));
                float3 diffuse = diff * _LightColor.rgb;
                
                float3 viewDir = normalize(-i.worldPos - _WorldSpaceCameraPos); 
                float3 reflectDir = reflect(-lightDir, normal);
                float spec = pow(max(dot(viewDir, reflectDir), 0.0), _SpecularPower);
                float3 specular = spec * _LightColor.rgb;
                
                float3 baseColor = tex2D(_MainTex, i.uv).rgb;
                
                float3 finalColor = ambient + diffuse + specular;
                finalColor *= baseColor;

                return float4(finalColor, 1.0);
            }
            ENDCG
        }
    }
    FallBack "Diffuse"
}
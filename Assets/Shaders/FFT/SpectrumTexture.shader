Shader "Ocean/HeightMapWithLighting"
{
    Properties
    {
        _MainTex("Height Map", 2D) = "white" {}  
        _Amplitude("Amplitude", Float) = 1.0     
        _Color("Base Color", Color) = (0.4, 0.6, 1, 1) 
    }

    SubShader
    {
        Tags { "RenderType" = "Opaque" }
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fwdbase

            #include <UnityLightingCommon.cginc>

            #include "UnityCG.cginc"

            sampler2D _MainTex;
            float _Amplitude;
            fixed4 _Color;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
            };

            struct v2f
            {
                float4 pos : SV_POSITION;
                float2 uv : TEXCOORD0;
                float3 worldNormal : NORMAL;
                float3 worldPos : TEXCOORD1;
            };
            
            v2f vert(appdata v)
            {
                v2f o;
                o.uv = v.uv;
                
                float height = tex2Dlod(_MainTex, float4(v.uv, 0, 0)).r * _Amplitude;
                float3 displacedVertex = v.vertex.xyz;
                displacedVertex.y += height; 

                o.pos = UnityObjectToClipPos(float4(displacedVertex, 1.0));
                o.worldPos = mul(unity_ObjectToWorld, float4(displacedVertex, 1.0)).xyz;
                o.worldNormal = UnityObjectToWorldNormal(v.normal);
                return o;
            }
            
            fixed4 frag(v2f i) : SV_Target
            {
                fixed3 worldNormal = normalize(i.worldNormal);
                fixed3 lightDir = normalize(_WorldSpaceLightPos0.xyz); 
                fixed3 viewDir = normalize(UnityWorldSpaceViewDir(i.worldPos));
                
                fixed3 diffuse = max(0, dot(worldNormal, lightDir)) * _LightColor0.rgb;
                fixed3 ambient = UNITY_LIGHTMODEL_AMBIENT.rgb;
                fixed3 halfDir = normalize(lightDir + viewDir);
                float specular = pow(max(0, dot(worldNormal, halfDir)), 16.0);

                fixed3 lighting = ambient + diffuse + specular;

                return fixed4(_Color.rgb * lighting, 1.0);
            }
            ENDCG
        }
    }
}
Shader "Custom/WaterbodySurfaceShader"
{
    Properties
    {
        _Color ("Color", Color) = (1,1,1,1)
        _MainTex ("Albedo (RGB)", 2D) = "white" {}
        _Glossiness ("Smoothness", Range(0,1)) = 0.5
        _Metallic ("Metallic", Range(0,1)) = 0.0
        _Alpha("Alpha", Range(0, 1)) = 1.0
        _HeightMap("Height Map", 2D) = "white" {}
        _FoamMap ("Foam Map", 2D) = "white" {}
        _FoamThreshold1("Foam Threshold 1", Float) = 0.5
        _FoamThreshold2("Foam Threshold 2", Float) = 0.5
        _FoamThreshold3("Foam Threshold 3", Float) = 0.5
        _Emission1("Emission 1", Float) = 0.0
        _Emission2("Emission 2", Float) = 0.0
        _Emission3("Emission 3", Float) = 0.0
        [ColorUsageAttribute(true, true)]_FoamColor1 ("Foam Color", Color) = (1,1,1,1)
        [ColorUsageAttribute(true, true)]_FoamColor2 ("Foam Color", Color) = (1,1,1,1)
        [ColorUsageAttribute(true, true)]_FoamColor3 ("Foam Color", Color) = (1,1,1,1)
        _DisplacementScale("Displacement Scale", Float) = 1.0
        _LightColor("Light Color", Color) = (1,1,1,1)
        _AmbientColor("Ambient Color", Color) = (0.2, 0.2, 0.2, 1) 
        _SpecularPower("Specular Power", Float) = 16.0
        _ChoppyFactor("Choppy Factor", Float) = 0.2
    }
    SubShader
    {
        Tags {
            "Queue"="Transparent"
            "RenderType"="Transparent"
        }
       // ZWrite Off
        
        Blend SrcAlpha OneMinusSrcAlpha
        LOD 200

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
        #pragma surface surf Standard alpha:premul fullforwardshadows vertex:vert

        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 3.0

        sampler2D _MainTex;

        struct Input
        {
            float2 uv_MainTex;
            float2 uv_HeightMap;
            float2 uv_FoamMap;
        };

        half _Glossiness;
        half _Metallic;
        fixed4 _Color;
        half _Alpha;

        sampler2D _HeightMap;
        sampler2D _FoamMap;
        float _DisplacementScale;
        float4 _LightColor;
        float4 _AmbientColor;
        float4 _FoamColor1, _FoamColor2, _FoamColor3;
        float _Emission1, _Emission2, _Emission3;
        float _ChoppyFactor;
        float _SpecularPower;
        float _FoamThreshold1, _FoamThreshold2, _FoamThreshold3;
        
        float4x4 _ObjectToWorld;
        float4x4 _WorldToObject;
        float4x4 _UnityMatrixVP;

        // Add instancing support for this shader. You need to check 'Enable Instancing' on materials that use the shader.
        // See https://docs.unity3d.com/Manual/GPUInstancing.html for more information about instancing.
        // #pragma instancing_options assumeuniformscaling
        UNITY_INSTANCING_BUFFER_START(Props)
            // put more per-instance properties here
        UNITY_INSTANCING_BUFFER_END(Props)

        void vert(inout appdata_full v)
        {
            float height = tex2Dlod(_HeightMap, float4(v.texcoord.xy, 0, 0)).g;
            float dx = tex2Dlod(_HeightMap, float4(v.texcoord.xy, 0, 0)).r;
            float dz = tex2Dlod(_HeightMap, float4(v.texcoord.xy, 0, 0)).b;

            v.vertex.y += height * _DisplacementScale;
            v.vertex.x += dx * _ChoppyFactor;
            v.vertex.z += dz * _ChoppyFactor;
        }
        
        void surf (Input IN, inout SurfaceOutputStandard o)
        {
            // Sample base texture
            float4 baseColor = _Color;

            // Calculate normal from height map
            float heightL = tex2D(_HeightMap, IN.uv_HeightMap + float2(-0.01, 0)).g;
            float heightR = tex2D(_HeightMap, IN.uv_HeightMap + float2(0.01, 0)).g;
            float heightD = tex2D(_HeightMap, IN.uv_HeightMap + float2(0, -0.01)).g;
            float heightU = tex2D(_HeightMap, IN.uv_HeightMap + float2(0, 0.01)).g;

            float foamValue = tex2D(_FoamMap, IN.uv_FoamMap).r;

            o.Albedo = _Color;
            // Foam Layer 1
            if (foamValue > _FoamThreshold1)
            {
                float foamIntensity1 = saturate((foamValue - _FoamThreshold1) / (1.0 - _FoamThreshold1));
                o.Albedo = lerp(o.Albedo, _FoamColor1.rgb, foamIntensity1);
                o.Emission += _FoamColor1.rgb * foamIntensity1 * _Emission1;
            }

            // Foam Layer 2
            if (foamValue > _FoamThreshold2)
            {
                float foamIntensity2 = saturate((foamValue - _FoamThreshold2) / (1.0 - _FoamThreshold2));
                o.Albedo = lerp(o.Albedo, _FoamColor2.rgb, foamIntensity2);
                o.Emission += _FoamColor2.rgb * foamIntensity2 * _Emission2;
            }

            // Foam Layer 3
            if (foamValue > _FoamThreshold3)
            {
                float foamIntensity3 = saturate((foamValue - _FoamThreshold3) / (1.0 - _FoamThreshold3));
                o.Albedo = lerp(o.Albedo, _FoamColor3.rgb, foamIntensity3);
                o.Emission += _FoamColor3.rgb * foamIntensity3 * _Emission3;
            }
            
            float3 normal = normalize(float3(heightL - heightR, 2.0, heightD - heightU));
            
            // Metallic and smoothness come from slider variables
            o.Metallic = _Metallic;
            o.Normal = normal;
            o.Smoothness = _Glossiness;
            o.Alpha = _Alpha;
        }
        ENDCG
    }
    FallBack "Diffuse"
}

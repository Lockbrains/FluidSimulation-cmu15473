
Shader "Custom/InstancedIndirectColor2D" {
    Properties {
        _ParticleRadius ("Particle Radius", Float) = 1.0 // Default value for particle radius
    }
    SubShader {
        Tags { "RenderType" = "Opaque" }

        Pass {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            struct appdata_t {
                float4 vertex   : POSITION;
                float4 color    : COLOR;
            };

            struct v2f {
                float4 vertex   : SV_POSITION;
                fixed4 color    : COLOR;
            }; 
            struct Particle
			{
                float pressure;
                float density;
                float2 currentForce;
                float2 velocity;
				float2 position;
                int cellId;
			};

            StructuredBuffer<Particle> Particles;
            float ParticleRadius;

            v2f vert(appdata_t i, uint instanceID: SV_InstanceID) {
                v2f o;

                // Scale the particle bc mesh is small
                float4x4 mat = {
                                    50*ParticleRadius, 0.f, 0.f, Particles[instanceID].position.x,
                                    0.f, 50*ParticleRadius, 0.f, Particles[instanceID].position.y,
                                    0.f, 0.f, 50*ParticleRadius, 0.f,
                                    0.f, 0.f, 0, 1.f
                                };

                float4 pos = mul(mat, i.vertex);
                o.vertex = UnityObjectToClipPos(pos);
                o.color = float4(0, 1, 1, 1);

                return o;
            }

            fixed4 frag(v2f i) : SV_Target {
                return i.color;
            }

            ENDCG
        }
    }
}

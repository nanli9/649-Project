//a shader for rendering
Shader "Custom/InstancedIndirectColor" {
	Properties
    {
        _Color ("Diffuse Color", Color) = (0.254, 0.419, 0.874, 1.0)
        _SpecColor ("Specular Color", Color) = (0.784, 0.863, 0.902, 1.0)
        _Shininess ("Shininess", Range(1, 100)) = 5.0
    }
    SubShader {
        Tags { "RenderType" = "Opaque" }
	
        Pass {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

			// Define properties
            fixed4 _Color;
            fixed4 _SpecColor;
            float _Shininess;

            struct appdata_t {
                float4 vertex   : POSITION;
                float4 color    : COLOR;
				float3 normal   : NORMAL;
            };

            struct v2f {
                float4 vertex   : SV_POSITION;
                fixed4 color    : COLOR;
				float3 normal   : NORMAL;
            }; 

            struct MeshProperties {
                float4x4 mat;
                float4 color;
            };

            StructuredBuffer<MeshProperties> _Properties;

            v2f vert(appdata_t i, uint instanceID: SV_InstanceID) {
                v2f o;

                float4 pos = mul(_Properties[instanceID].mat, i.vertex);
                o.vertex = UnityObjectToClipPos(pos);

				o.normal = UnityObjectToWorldNormal(i.normal);

                o.color = _Properties[instanceID].color;

                return o;
            }

            fixed4 frag(v2f i) : SV_Target {
				float3 N = normalize(i.normal);
                float3 L = normalize(_WorldSpaceLightPos0.xyz);
                float3 V = normalize(_WorldSpaceCameraPos - i.vertex);

				float NdotL = max(0, dot(N, L));
                float3 diffuse = _Color.rgb * NdotL;

                float3 R = reflect(-L, N);
                float spec = pow(max(dot(V, R), 0.0), _Shininess);
                float3 specular = _SpecColor.rgb * spec;

                // Combine lighting components
                return fixed4(diffuse + specular, 1.0);
                //return i.color;
            }

            ENDCG
        }
    }
}
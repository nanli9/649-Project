Shader "Unlit/RayMarching"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // make fog work
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }
			Texture3D<float4> DensityMap;
			SamplerState samplerDensityMap; // Define a sampler state
			const float3 boundsSize;
			static const float TinyNudge = 0.01;
			const float densityMultiplier;
			const float stepSize;
			const float volumeValueOffset;

			float2 RayBox(float3 boundsMin, float3 boundsMax, float3 rayOrigin, float3 rayDir)
			{
				float3 invRayDir = 1 / rayDir;

                float3 t0 = (boundsMin - rayOrigin) * invRayDir;
                float3 t1 = (boundsMax - rayOrigin) * invRayDir;
                float3 tmin = min(t0, t1);
                float3 tmax = max(t0, t1);

                float dstA = max(max(tmin.x, tmin.y), tmin.z);
                float dstB = min(tmax.x, min(tmax.y, tmax.z));

                // CASE 1: ray intersects box from outside (0 <= dstA <= dstB)
                // dstA is dst to nearest intersection, dstB dst to far intersection

                // CASE 2: ray intersects box from inside (dstA < 0 < dstB)
                // dstA is the dst to intersection behind the ray, dstB is dst to forward intersection

                // CASE 3: ray misses box (dstA > dstB)

                float dstToBox = max(0, dstA);
                float dstInsideBox = max(0, dstB - dstToBox);
                return float2(dstToBox, dstInsideBox);
			}

			float SampleDensity(float3 pos)
            {
                float3 uvw = pos  / boundsSize;

                const float epsilon = 0.0001;
                bool isEdge = any(uvw >= 1 - epsilon || uvw <= epsilon);
                if (isEdge) return -volumeValueOffset;

                return DensityMap.SampleLevel(samplerDensityMap, uvw, 0).r - volumeValueOffset;
            }

			float3 RayMarch(float2 uv, float stepSize)
			{
				float3 localViewVector = mul(unity_CameraInvProjection, float4(uv * 2 - 1, 0, -1));
				float3 rayDir = normalize(mul(unity_CameraToWorld, float4(localViewVector, 0)));

				float2 boundsDstInfo = RayBox(0, boundsSize, _WorldSpaceCameraPos, rayDir);
				float3 entryPoint = _WorldSpaceCameraPos + rayDir * (boundsDstInfo[0] + TinyNudge);

				float densityAlongViewRay = 0;

				for (float dst = 0; dst < boundsDstInfo[1] - TinyNudge * 2; dst += stepSize)
				{
					float3 samplePos = entryPoint + rayDir * dst;
					float densityAlongStep = SampleDensity(samplePos) + densityMultiplier * stepSize;
					densityAlongViewRay += densityAlongStep;
				}


				return densityAlongViewRay;
			}
            fixed4 frag (v2f i) : SV_Target
            {

                return float4(RayMarch(i.uv, stepSize), 1.0);
            }
            ENDCG
        }
    }
}

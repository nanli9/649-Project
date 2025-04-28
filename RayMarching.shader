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
			const float3 dirToSun;
			const float3 fluidColor;
			const float extinctionCoeff;
			static const float TinyNudge = 0.01;
			const float densityMultiplier;
			const float stepSize;
			const float volumeValueOffset;
			const float indexOfRefraction = 1.333f;
			const float iorAir = 1.0f;
			float3 hapticPointPos;
            float hapticPointRadius;

			float2 RayBox(float3 boundsMin, float3 boundsMax, float3 rayOrigin, float3 rayDir)
			{
				float3 invRayDir = 1 / rayDir;

				//float offset = 1;
				//boundsMin += float3(offset,offset,offset);
				//boundsMax -= float3(offset,offset,offset);

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
			float RaySphereIntersect(float3 sphereCenter, float sphereRadius, float3 rayOrigin, float3 rayDir) {
                float3 oc = rayOrigin - sphereCenter;

                float b = 2.0 * dot(oc, rayDir);
                float c = dot(oc, oc) - sphereRadius * sphereRadius;
                float discriminant = b * b - 4 * c; // 4*a*c -> 4*c

                if (discriminant < 0) {
                    return -1.0; // No intersection
                }
                else {
                    float sqrtDiscriminant = sqrt(discriminant);
                    float t1 = (-b - sqrtDiscriminant) / 2.0;
                    float t2 = (-b + sqrtDiscriminant) / 2.0;
                    if (t1 >= 0) return t1; // Closest intersection in front
                    if (t2 >= 0) return t2; // Ray origin is inside sphere
                    return -1.0; // Both intersections are behind
                }
            }
			float3 Light(float3 pos, float3 dir)
            {
                //return SampleEnvironmentAA(pos, dir);
				const float3 colGround = float3(0.35, 0.3, 0.35) * 0.53;
                const float3 colSkyHorizon = float3(1, 1, 1);
                const float3 colSkyZenith = float3(0.08, 0.37, 0.73);

                float sun = pow(max(0, dot(dir, dirToSun)), 500) * 1;
                float skyGradientT = pow(smoothstep(0, 0.4, dir.y), 0.35);
                float groundToSkyT = smoothstep(-0.01, 0, dir.y);
                float3 skyGradient = lerp(colSkyHorizon, colSkyZenith, skyGradientT);

                return lerp(colGround, skyGradient, groundToSkyT) + sun * (groundToSkyT >= 1);
            }
			float SampleDensity(float3 pos)
            {
                float3 uvw = pos  / boundsSize;

                const float epsilon = 0.0001;
                bool isEdge = any(uvw >= 1 - epsilon || uvw <= epsilon);
                if (isEdge) return -volumeValueOffset;

                return DensityMap.SampleLevel(samplerDensityMap, uvw, 0).r - volumeValueOffset;
            }
			float CalculateDensityAlongRay(float3 rayPos, float3 rayDir, float stepSize)
            {
                // Test for non-normalize ray and return 0 in that case.
                // This happens when refract direction is calculated, but ray is totally reflected
                if (dot(rayDir, rayDir) < 0.9) return 0;

                float2 boundsDstInfo = RayBox(float3(0,0,0), boundsSize, rayPos, rayDir);
                float dstToBounds = boundsDstInfo[0];
                float dstThroughBounds = boundsDstInfo[1];
                if (dstThroughBounds <= 0) return 0;

                float dstTravelled = 0;
                float opticalDepth = 0;
                float nudge = stepSize * 0.5;
                float3 entryPoint = rayPos + rayDir * (dstToBounds + nudge);
                dstThroughBounds -= (nudge + TinyNudge);

                while (dstTravelled < dstThroughBounds)
                {
                    rayPos = entryPoint + rayDir * dstTravelled;
                    float density = SampleDensity(rayPos) * densityMultiplier * stepSize;
                    if (density > 0)
                    {
                        opticalDepth += density;
                    }
                    dstTravelled += stepSize;
                }

                return opticalDepth;
            }
			uint NextRandom(inout uint state)
            {
                state = state * 747796405 + 2891336453;
                uint result = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
                result = (result >> 22) ^ result;
                return result;
            }
			float RandomValue(inout uint state)
            {
                return NextRandom(state) / 4294967295.0; // 2^32 - 1
            }
			float3 Transmittance(float thickness)
            {
                return exp(-thickness * extinctionCoeff);
            }
			struct LightResponse
            {
                float3 reflectDir;
                float3 refractDir;
                float reflectWeight;
                float refractWeight;
            };
			struct SurfaceInfo
            {
                float3 pos;
                float3 normal;
                float densityAlongRay;
                bool foundSurface;
				bool isHapticPoint;
            };
			bool IsInsideFluid(float3 pos)
            {
                // Check if the position is within the bounds of the fluid volume (starting at origin)
				bool withinBounds = all(pos >= float3(0, 0, 0)) && all(pos <= boundsSize);

				// Also check if there's density at that position
				return withinBounds && SampleDensity(pos) > 0;
            }
			float3 Refract(float3 inDir, float3 normal, float iorA, float iorB)
            {
                float refractRatio = iorA / iorB;
                float cosAngleIn = -dot(inDir, normal);
                float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
                if (sinSqrAngleOfRefraction > 1) return 0; // Ray is fully reflected, no refraction occurs

                float3 refractDir = refractRatio * inDir + (refractRatio * cosAngleIn - sqrt(1 - sinSqrAngleOfRefraction)) * normal;
                return refractDir;
            }
			float3 Reflect(float3 inDir, float3 normal)
            {
                return inDir - 2 * dot(inDir, normal) * normal;
            }
			float3 CalculateClosestFaceNormal(float3 boxSize, float3 p)
            {
                float3 oMin = abs(p - float3(0, 0, 0)); // Distance to the faces at the origin
				float3 oMax = abs(p - boxSize);         // Distance to the faces at boxSize

				if (oMin.x < oMax.x && oMin.x < oMin.y && oMin.x < oMax.y && oMin.x < oMin.z && oMin.x < oMax.z)
					return float3(-1, 0, 0); // Closest to the x=0 face
				if (oMax.x < oMin.y && oMax.x < oMax.y && oMax.x < oMin.z && oMax.x < oMax.z)
					return float3(1, 0, 0);  // Closest to the x=boxSize.x face
				if (oMin.y < oMax.y && oMin.y < oMin.z && oMin.y < oMax.z)
					return float3(0, -1, 0); // Closest to the y=0 face
				if (oMax.y < oMin.z && oMax.y < oMax.z)
					return float3(0, 1, 0);  // Closest to the y=boxSize.y face
				if (oMin.z < oMax.z)
					return float3(0, 0, -1); // Closest to the z=0 face
				return float3(0, 0, 1);  // Closest to the z=boxSize.z face

            }
			float3 CalculateNormal(float3 pos)
            {
                float3 uvw = pos / boundsSize;

                const float s = 0.1;
                float3 offsetX = float3(1, 0, 0) * s;
                float3 offsetY = float3(0, 1, 0) * s;
                float3 offsetZ = float3(0, 0, 1) * s;

                float dx = SampleDensity(pos - offsetX) - SampleDensity(pos + offsetX);
                float dy = SampleDensity(pos - offsetY) - SampleDensity(pos + offsetY);
                float dz = SampleDensity(pos - offsetZ) - SampleDensity(pos + offsetZ);

                float3 volumeNormal = normalize(float3(dx, dy, dz));

                // Smoothly flatten normals out at boundary edges
                float3 o = boundsSize - pos;
                float faceWeight = min(o.x, min(o.y, o.z));
                float3 faceNormal = CalculateClosestFaceNormal(boundsSize, pos);
                const float smoothDst = 0.3;
                const float smoothPow = 5;
                faceWeight = (1 - smoothstep(0, smoothDst, faceWeight)) * (1 - pow(saturate(volumeNormal.y), smoothPow));

                return normalize(volumeNormal * (1 - faceWeight) + faceNormal * (faceWeight));
            }
			float CalculateReflectance(float3 inDir, float3 normal, float iorA, float iorB)
            {
                float refractRatio = iorA / iorB;
                float cosAngleIn = -dot(inDir, normal);
                float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
                if (sinSqrAngleOfRefraction >= 1) return 1; // Ray is fully reflected, no refraction occurs

                float cosAngleOfRefraction = sqrt(1 - sinSqrAngleOfRefraction);
                // Perpendicular polarization
                float rPerpendicular = (iorA * cosAngleIn - iorB * cosAngleOfRefraction) / (iorA * cosAngleIn + iorB * cosAngleOfRefraction);
                rPerpendicular *= rPerpendicular;
                // Parallel polarization
                float rParallel = (iorB * cosAngleIn - iorA * cosAngleOfRefraction) / (iorB * cosAngleIn + iorA * cosAngleOfRefraction);
                rParallel *= rParallel;

                // Return the average of the perpendicular and parallel polarizations
                return (rPerpendicular + rParallel) / 2;
            }
			LightResponse CalculateReflectionAndRefraction(float3 inDir, float3 normal, float iorA, float iorB)
            {
                LightResponse result;

                result.reflectWeight = CalculateReflectance(inDir, normal, iorA, iorB);
                result.refractWeight = 1 - result.reflectWeight;

                result.reflectDir = Reflect(inDir, normal);
                result.refractDir = Refract(inDir, normal, iorA, iorB);

                return result;
            }
			SurfaceInfo FindNextSurface(float3 origin, float3 rayDir, bool findNextFluidEntryPoint, uint rngState, float rngWeight, float maxDst)
            {
                SurfaceInfo info = (SurfaceInfo)0;
				info.isHapticPoint = false;

                if (dot(rayDir, rayDir) < 0.5) return info;

				float sphereHitDistance = RaySphereIntersect(hapticPointPos, hapticPointRadius, origin, rayDir); //returns -1 if no hit.

				// First check sphere intersection
				if (sphereHitDistance > 0 && sphereHitDistance < maxDst) { //hit
					info.pos = origin + rayDir * sphereHitDistance;
					info.foundSurface = true;
					info.densityAlongRay = 0; //reset density
					info.isHapticPoint = true;
					return info;
				}

                float2 boundsDstInfo = RayBox(float3(0,0,0), boundsSize, origin, rayDir);
                float r = (RandomValue(rngState) - 0.5) * stepSize * 0.4 * 1;
                bool hasExittedFluid = !IsInsideFluid(origin);
                origin = origin + rayDir * (boundsDstInfo.x + r);

                //float stepSize = viewMarchStepSize;
                bool hasEnteredFluid = false;
                float3 lastPosInFluid = origin;

                float dstToTest = boundsDstInfo[1] - (TinyNudge) * 2;

                for (float dst = 0; dst < dstToTest; dst += stepSize)
                {
                    bool isLastStep = dst + stepSize >= dstToTest;
                    float3 samplePos = origin + rayDir * dst;
                    float thickness = SampleDensity(samplePos) * densityMultiplier * stepSize;
                    bool insideFluid = thickness > 0;
                    if (insideFluid)
                    {
                        hasEnteredFluid = true;
                        lastPosInFluid = samplePos;
                        if (dst <= maxDst)
                        {
                            info.densityAlongRay += thickness;
                        }
                    }

                    if (!insideFluid) hasExittedFluid = true;

                    bool found;
                    if (findNextFluidEntryPoint) found = insideFluid && hasExittedFluid;
                    else found = hasEnteredFluid && (!insideFluid || isLastStep);

                    if (found)
                    {
                        info.pos = lastPosInFluid;
                        info.foundSurface = true;
                        break;
                    }
                }

                return info;
            }

			float3 RayMarch(float2 uv, float stepSize)
			{
				uint rngState = (uint)(uv.x * 1243 + uv.y * 96456);
				float3 localViewVector = mul(unity_CameraInvProjection, float4(uv * 2 - 1, 0, -1));
				float3 rayDir = normalize(mul(unity_CameraToWorld, float4(localViewVector, 0)));
				float3 rayPos = _WorldSpaceCameraPos.xyz;
                bool travellingThroughFluid = IsInsideFluid(rayPos);

				float3 transmittance = 1;
                float3 light = 0;

				int numRefractions = 2;
				float lightStepSize = 0.02;

				float3 sphereColor = float3(1,1,1);
				for (int i = 0; i < numRefractions; i++)
				{
					float densityStepSize = lightStepSize * (i + 1); // increase step size with each iteration
					bool searchForNextFluidEntryPoint = !travellingThroughFluid;

					SurfaceInfo surfaceInfo = FindNextSurface(rayPos, rayDir, searchForNextFluidEntryPoint, rngState, i == 0 ? 1 : 0, 1e10); // Removed cubeHit.dst

					if (!surfaceInfo.foundSurface) break;
					if (surfaceInfo.isHapticPoint)
						return sphereColor; // Use sphere color
					transmittance *= Transmittance(surfaceInfo.densityAlongRay) * fluidColor;

					float3 normal = CalculateNormal(surfaceInfo.pos);
					if (dot(normal, rayDir) > 0) normal = -normal;

					//return normal;

					// Indicies of refraction
					float iorA = travellingThroughFluid ? indexOfRefraction : iorAir;
					float iorB =travellingThroughFluid ? iorAir : indexOfRefraction;

					// Calculate reflection and refraction, and choose which path to follow
					LightResponse lightResponse = CalculateReflectionAndRefraction(rayDir, normal, iorA, iorB);
					float densityAlongRefractRay = CalculateDensityAlongRay(surfaceInfo.pos, lightResponse.refractDir, densityStepSize);
					float densityAlongReflectRay = CalculateDensityAlongRay(surfaceInfo.pos, lightResponse.reflectDir, densityStepSize);
					bool traceRefractedRay = densityAlongRefractRay * lightResponse.refractWeight > densityAlongReflectRay * lightResponse.reflectWeight;
					travellingThroughFluid = traceRefractedRay != travellingThroughFluid;

					// Approximate less interesting path
					if (traceRefractedRay) light += Light(surfaceInfo.pos, lightResponse.reflectDir) * transmittance * Transmittance(densityAlongReflectRay) * lightResponse.reflectWeight;
					else light += Light(surfaceInfo.pos, lightResponse.refractDir) * transmittance * Transmittance(densityAlongRefractRay) * lightResponse.refractWeight;

					// Set up ray for more interesting path
					rayPos = surfaceInfo.pos;
					rayDir = traceRefractedRay ? lightResponse.refractDir : lightResponse.reflectDir;
					transmittance *= (traceRefractedRay ? lightResponse.refractWeight : lightResponse.reflectWeight);
				}

				// Approximate remaining path
				float densityRemainder = CalculateDensityAlongRay(rayPos, rayDir, lightStepSize);
				light += Light(rayPos, rayDir) * transmittance * Transmittance(densityRemainder);

				return light;
			}

            fixed4 frag (v2f i) : SV_Target
            {

                return float4(RayMarch(i.uv, stepSize), 1.0);
            }
            ENDCG
        }
    }
}

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static UnityEngine.GraphicsBuffer;

public class FluidRenderer : MonoBehaviour
{
    [Header("Setting")]
    public float densityOffset;
    public float stepSize;
    public float densityMultiplier;
    public float extinctionCoeff;
    public Vector3 dirToSun;
    public Vector3 fluidColor;

    [Header("Reference")]
    public Shader shader;
    public PBF pbf;

    Material raymarchMat;

    void Start()
    {
        raymarchMat = new Material(shader);
        Camera.main.depthTextureMode = DepthTextureMode.Depth;
    }

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        if (pbf.DensityMap != null)
        {
            SetShaderParams();
            Graphics.Blit(source, destination, raymarchMat);
        }
        else
        {
            Graphics.Blit(source, destination);
        }
    }

    void SetShaderParams()
    {
        raymarchMat.SetTexture("DensityMap", pbf.DensityMap);
        float boundsX = pbf.boundingBox_x.y - pbf.boundingBox_x.x;
        float boundsY = pbf.boundingBox_y.y - pbf.boundingBox_y.x;
        float boundsZ = pbf.boundingBox_z.y - pbf.boundingBox_z.x;

        Vector3 BoundsSize = new Vector3(boundsX, boundsY, boundsZ);
        raymarchMat.SetVector("boundsSize", BoundsSize);
        raymarchMat.SetVector("dirToSun", dirToSun);
        raymarchMat.SetVector("fluidColor", fluidColor);
        raymarchMat.SetFloat("stepSize", stepSize);
        raymarchMat.SetFloat("densityMultiplier", densityMultiplier);
        raymarchMat.SetFloat("volumeValueOffset", densityOffset);
        raymarchMat.SetFloat("extinctionCoeff", extinctionCoeff);

        raymarchMat.SetFloat("hapticPointRadius", 0.5f);
        //raymarchMat.SetVector("hapticPointPos", pbf.hapticInteractionPoint);
        raymarchMat.SetVector("hapticPointPos", new Vector3(5,10,5));
    }
}

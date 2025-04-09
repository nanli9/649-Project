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
        raymarchMat.SetFloat("stepSize", stepSize);
        raymarchMat.SetFloat("densityMultiplier", densityMultiplier);
        raymarchMat.SetFloat("volumeValueOffset", densityOffset);
    }
}

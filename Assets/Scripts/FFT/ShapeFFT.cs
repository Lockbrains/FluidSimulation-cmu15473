using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using Random = UnityEngine.Random;
using Vector2 = UnityEngine.Vector2;

/// <summary>
/// For generating the Wave Shape using FFT.
/// </summary>
public class ShapeFFT : MonoBehaviour
{
    [SerializeField] private int resolution = 128;
    [SerializeField] private float frequency = 5f;
    [SerializeField] private float amplitude = 1f;
    [SerializeField] private float speed = 1.0f;

    [Header("For Spectrum Visualization")]
    [SerializeField] private MeshRenderer spectrumPlaneMR;
    
    [Header("Phillips Spectrum")]
    [SerializeField] private float windSpeed = 10f; // wind speed
    [SerializeField] private float g = 9.81f;       // gravity accelaration
    [SerializeField] private float A = 0.01f;       // adjustment factor
     
    private float[,] spectrum;
    private Complex[,] fftOutput;
    private int spectrumSize;
    private Texture2D spectrumTexture;
    
    [SerializeField] private ComputeShader computeShader;
    private int kernelHandle;
    private int kernelInitSpectrum;
    private int kernelUpdatePhase;
    private int kernelFFT1D_Horizontal;
    private int kernelFFT1D_Vertical;
    private int kernelIFFT1D_Horizontal;
    private int kernelIFFT1D_Vertical;
    private int kernelFinalOutput;

    private RenderTexture tempBuffer;
    private RenderTexture heightMap;
    
    private int threadGroupX, threadGroupY;

    void Awake()
    {
        threadGroupX = Mathf.CeilToInt((float)resolution / 8);
        threadGroupY = Mathf.CeilToInt((float)resolution / 8);
        
        kernelInitSpectrum = computeShader.FindKernel("InitSpectrum");
        kernelUpdatePhase = computeShader.FindKernel("UpdatePhase");
        kernelFFT1D_Horizontal = computeShader.FindKernel("FFT1D_Horizontal");
        kernelFFT1D_Vertical = computeShader.FindKernel("FFT1D_Vertical");
        kernelIFFT1D_Horizontal = computeShader.FindKernel("IFFT1D_Horizontal");
        kernelIFFT1D_Vertical = computeShader.FindKernel("IFFT1D_Vertical");
        kernelFinalOutput = computeShader.FindKernel("FinalOutput");
    }
    
    // Start is called before the first frame update
    void Start()
    {
        tempBuffer = new RenderTexture(resolution, resolution, 0, RenderTextureFormat.ARGBFloat)
        {
            enableRandomWrite = true
        };
        tempBuffer.Create();

        heightMap = new RenderTexture(resolution, resolution, 0, RenderTextureFormat.RFloat)
        {
            enableRandomWrite = true
        };
        heightMap.Create();

        computeShader.SetInt("_Resolution", resolution);
        computeShader.SetTexture(kernelInitSpectrum, "_TempBuffer", tempBuffer);
        computeShader.Dispatch(kernelInitSpectrum, threadGroupX, threadGroupY, 1);
        
        // spectrum = new float[resolution, resolution];
        // spectrumTexture = new Texture2D(resolution, resolution);
        
        //InitializeSpectrum();
    }

    // Update is called once per frame
    void Update()
    {
        computeShader.SetInt("_Resolution", resolution);
        computeShader.SetFloat("_Time", Time.time * speed);
        
        computeShader.SetTexture(kernelInitSpectrum, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelUpdatePhase, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelFFT1D_Horizontal, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelFFT1D_Vertical, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelIFFT1D_Horizontal, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelIFFT1D_Vertical, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelFinalOutput, "_TempBuffer", tempBuffer);
        computeShader.SetTexture(kernelFinalOutput, "_HeightMap", heightMap);


        computeShader.Dispatch(kernelUpdatePhase, threadGroupX, threadGroupY, 1);
        computeShader.Dispatch(kernelFFT1D_Horizontal, threadGroupX, threadGroupY, 1);
        computeShader.Dispatch(kernelFFT1D_Vertical, threadGroupX, threadGroupY, 1);
        computeShader.Dispatch(kernelIFFT1D_Horizontal, threadGroupX, threadGroupY, 1);
        computeShader.Dispatch(kernelIFFT1D_Vertical, threadGroupX, threadGroupY, 1);
        
        computeShader.Dispatch(kernelFinalOutput, threadGroupX, threadGroupY, 1);
        
        spectrumPlaneMR.material.SetTexture("_MainTex", heightMap);
        
        RenderTexture.active = heightMap;
        
        Texture2D tempTex = new Texture2D(resolution, resolution, TextureFormat.RFloat, false);
        tempTex.ReadPixels(new Rect(0, 0, resolution, resolution), 0, 0);
        tempTex.Apply();
        
        float debugValue = tempTex.GetPixel(0, 0).r; 
        Debug.Log($"Pixel (0,0) Value: {debugValue}");
    }
    
    private void InitializeSpectrum()
    {
        for (int x = 0; x < resolution; x++)
        {
            for (int y = 0; y < resolution; y++)
            {
                // wave vectors
                float kx = (x - resolution / 2) * 2.0f * Mathf.PI / resolution;
                float kz = (y - resolution / 2) * 2.0f * Mathf.PI / resolution;
                float k = Mathf.Sqrt(kx * kx + kz * kz);

                if (k == 0) continue;
                
                // wind direction
                Vector2 windDir = new Vector2(1, 0);
                float windFactor = Vector2.Dot(new Vector2(kx, kz).normalized, windDir);
                
                // Phillips Spectrum
                float L = windSpeed * windSpeed / g;
                float damping = Mathf.Exp(-1.0f / (k * L * k * L));
                float phillips = A * damping / (k * k * k * k) * windFactor * windFactor; 
                
                // Time Factor
                float omega = Mathf.Sqrt(g * k);
                float phase = omega * Time.time;
                float realPart = Mathf.Sqrt(phillips) * Mathf.Cos(phase);
                
                spectrum[x, y] = realPart; 
            }
        }

        GenerateTextureFromSpectrum();
    }

    private void GenerateTextureFromSpectrum()
    {
        int size = spectrum.GetLength(0);

        if (!spectrumTexture)
        {
            spectrumTexture = new Texture2D(size, size);
        }
        
        Color[] pixels = new Color[size * size];
        
        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                float value = spectrum[x, y];
                pixels[x + y * size] = new Color(value, value, value, 1.0f);
            }
        }
        
        spectrumTexture.SetPixels(pixels);
        spectrumTexture.Apply();
    }

    public void GenerateSpectrum()
    {
        spectrumPlaneMR.material.mainTexture = spectrumTexture;
    }

    private void PerformFFT()
    {
        fftOutput = FFT2D(spectrum);
        
        // int size = fftOutput.GetLength(0);
        //
        // for (int x = 0; x < size; x++)
        // {
        //     for (int y = 0; y < size; y++)
        //     {
        //         float magnitude = (float) fftOutput[x, y].Magnitude;
        //         magnitude = Mathf.Clamp01(magnitude / 10f);
        //         
        //         spectrumTexture.SetPixel(x, y, new Color(magnitude, magnitude, magnitude, 1.0f));
        //     }
        // }
        //
        // spectrumTexture.Apply();
    }
    
    private void PerformIFFT()
    {
        Complex[,] ifftOutput = IFFT2D(fftOutput);
        
        int size = ifftOutput.GetLength(0);
        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                float height = (float)ifftOutput[x, y].Real; 
                height = Mathf.Clamp01((height + 1f) / 2f);
                
                spectrumTexture.SetPixel(x, y, new Color(height, height, height, 1.0f));
            }
        }
        
        spectrumTexture.Apply();
    }
    
    #region 1D FFT
    public static Complex[] FFT1D(Complex[] input)
    {
        int N = input.Length;
        if (N <= 1) return input;

        var even = new Complex[N / 2];
        var odd  = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = input[i * 2];
            odd [i] = input[i * 2 + 1];
        }

        var fftEven = FFT1D(even);
        var fftOdd = FFT1D(odd);

        var result = new Complex[N];
        for (int k = 0; k < N / 2; k++)
        {
            Complex t = Complex.FromPolarCoordinates(1.0, -2.0 * Math.PI * k / N) * fftOdd[k];
            result[k] = fftEven[k] + t;
            result[k + N / 2] = fftEven[k] - t;
        }

        return result;
    }

    #endregion
    
    #region 2D FFT

    private Complex[,] FFT2D(float[,] input)
    {
        int size = input.GetLength(0);
        Complex[,] complexData = new Complex[size, size];

        // initialize input data
        for (int x = 0; x < size; x++)
        {
            for (int y = 0; y < size; y++)
            {
                complexData[x, y] = new Complex(input[x, y], 0);
            }
        }
        
        // run FFT1D for every row
        for (int y = 0; y < size; y++)
        {
            Complex[] row = new Complex[size];
            for (int x = 0; x < size; x++)
            {
                row[x] = complexData[x, y];
            }
            
            Complex[] fftRow = FFT1D(row);
            for (int x = 0; x < size; x++)
            {
                complexData[x, y] = fftRow[x];
            }
        }
        
        // run FFT1D for every col
        for (int x = 0; x < size; x++)
        {
            Complex[] column = new Complex[size];
            for (int y = 0; y < size; y++)
            {
                column[y] = complexData[x, y];
            }
            
            Complex[] fftColumn = FFT1D(column);
            for (int y = 0; y < size; y++)
            {
                complexData[x, y] = fftColumn[y]; 
            }
        }
        
        return complexData;
    }
    
    #endregion

    #region Inverse FFT 1D
    public static Complex[] IFFT1D(Complex[] input)
    {
        int N = input.Length;
        if (N <= 1) return input;

        var even = new Complex[N / 2];
        var odd = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++)
        {
            even[i] = input[i * 2];
            odd[i] = input[i * 2 + 1];
        }

        var ifftEven = IFFT1D(even);
        var ifftOdd = IFFT1D(odd);

        var result = new Complex[N];
        for (int k = 0; k < N / 2; k++)
        {
            Complex t = Complex.FromPolarCoordinates(1.0, +2.0 * Math.PI * k / N) * ifftOdd[k];
            result[k] = (ifftEven[k] + t) / 2.0;
            result[k + N / 2] = (ifftEven[k] - t) / 2.0;
        }

        return result;
    }

    #endregion
    
    #region Inverse FFT 2D
    private Complex[,] IFFT2D(Complex[,] input)
    {
        int size = input.GetLength(0);
        Complex[,] complexData = new Complex[size, size];

        // Inverse 1D FFT each row
        for (int y = 0; y < size; y++)
        {
            Complex[] row = new Complex[size];
            for (int x = 0; x < size; x++)
            {
                row[x] = input[x, y];
            }

            Complex[] ifftRow = IFFT1D(row);
            for (int x = 0; x < size; x++)
            {
                complexData[x, y] = ifftRow[x];
            }
        }

        // Inverse 1D FFT each column
        for (int x = 0; x < size; x++)
        {
            Complex[] column = new Complex[size];
            for (int y = 0; y < size; y++)
            {
                column[y] = complexData[x, y];
            }

            Complex[] ifftColumn = IFFT1D(column);
            for (int y = 0; y < size; y++)
            {
                complexData[x, y] = ifftColumn[y];
            }
        }

        return complexData;
    }
    #endregion
} 

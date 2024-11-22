using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using Unity.VisualScripting;
using UnityEngine;
using Vector2 = UnityEngine.Vector2;

/// <summary>
/// For generating the Wave Shape using FFT.
/// Creator: Lingheng Tony Tao
/// </summary>
public class ShapeFFT : MonoBehaviour
{
    [SerializeField] private int N = 256; // N
    [SerializeField] private float L = 1000f; // Patch Size
    [SerializeField] private float A = 0.01f; // Adjustment
    [SerializeField] private float V = 40f; // Wind Speed
    [SerializeField] private Vector2 windDirection = new Vector2(1, 1).normalized; 
    [SerializeField] private float frequency = 5f;

    [Header("For Spectrum Visualization")]
    [SerializeField] private MeshRenderer h0kRenderer;
    [SerializeField] private MeshRenderer h0MinusKRenderer;
    [SerializeField] private MeshRenderer butterflyTextureRenderer;
    [SerializeField] private MeshRenderer updateRenderer;
    [SerializeField] private MeshRenderer FFT1DResultTexture;
    
    // Private Variables 
    private readonly float gravity = 9.81f;

    private Complex[,] h0k;
    private Complex[,] h0MinusK;

    private Texture2D h0kTexture;
    private Texture2D h0MinusKTexture;
    private Texture2D butterflyTexture;

    private float time = 0f;

    void Start()
    {
        InitializeSpectrum();

        int[] bitReversedIndices = GenerateBitReversedIndices(N);
        InitializeButterflyTexture(N, bitReversedIndices);
        
        butterflyTextureRenderer.material.SetTexture("_MainTex", butterflyTexture);
        butterflyTexture.filterMode = FilterMode.Point;
        //butterflyTexture.wrapMode = TextureWrapMode.Repeat;
        //butterflyTextureRenderer.material.SetTextureScale("_MainTex", new Vector2(1.0f, 1.0f / (Mathf.Log(N, 2))));
        
        h0kTexture = GenerateTexture(h0k);
        h0MinusKTexture = GenerateTexture(h0MinusK);
        
        h0kRenderer.material.SetTexture("_MainTex", h0kTexture);
        h0MinusKRenderer.material.SetTexture("_MainTex", h0MinusKTexture);
    }
    // Update is called once per frame
    void Update()
    {
        // Update Fourier Components
        time += Time.deltaTime;

        Texture2D hktTexture = UpdateFourierComponents(h0k, h0MinusK, time);
        
        //updateRenderer.material.SetTexture("_MainTex", hktTexture); 
        
        Complex[,] hktData = ExtractDataFromTexture(hktTexture);
        
        // Visualize 1D FFT Result
        Texture2D fftResultTexture = GenerateTexture(hktData);
        FFT1DResultTexture.material.SetTexture("_MainTex", fftResultTexture);
    }
    
    #region Phase 1 Initialization
    // Phase 1: Produces the spectrum textures containing \Tilde{h}_0(\bold{k}) and \Tilde{h}_0^*(-\bold{k})
    // Page 24, 4.2.3
    private (float, float) GenerateGaussianRandom()
    {
        float u1 = UnityEngine.Random.value;
        float u2 = UnityEngine.Random.value;

        float z0 = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Cos(2.0f * Mathf.PI * u2);
        float z1 = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Sin(2.0f * Mathf.PI * u2);

        return (z0, z1);
    }
    
    
    private void InitializeSpectrum()
    {
        float windLength = V * V / gravity; // L = V^2/g;
        h0k = new Complex[N, N];
        h0MinusK = new Complex[N,N];
        
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Wave Vector \bold{k}
                // Equation (4.2)
                float kx = (2 * Mathf.PI * i - Mathf.PI * N ) / L;
                float kz = (2 * Mathf.PI * j - Mathf.PI * N ) / L;
                Vector2 k = new Vector2(kx, kz);
                
                // Clamp the size 
                float kLength = k.magnitude;
                if (kLength < 1e-6f) kLength = 1e-6f;

                // Phillips Spectrum
                // Equation (3.5)
                float kDotW = Vector2.Dot(k.normalized, windDirection);
                float kDotW2 = kDotW * kDotW;
                float phillips = A / Mathf.Pow(kLength, 4) * Mathf.Exp(-1 / Mathf.Pow(kLength * windLength, 2));
                phillips *= kDotW2;

                // Suppressing very small waves
                if (kLength < 1e-3f)
                    phillips = 0;

                // Generate Gaussian Random Numbers
                var (gaussianReal, gaussianImag) = GenerateGaussianRandom();

                // Initialize h0(k) and h0(-k)
                h0k[i, j] = new Complex(gaussianReal * Mathf.Sqrt(phillips / 2), gaussianImag * Mathf.Sqrt(phillips / 2));
                h0MinusK[i, j] = new Complex(h0k[i, j].Real, -h0k[i, j].Imaginary);
            }
        }
    }

    int[] GenerateBitReversedIndices(int N)
    {
        int[] indices = new int[N];
        int log2N = Mathf.RoundToInt(Mathf.Log(N, 2));

        for (int i = 0; i < N; i++)
        {
            int reversed = 0;
            int index = i;

            for (int j = 0; j < log2N; j++)
            {
                reversed = (reversed << 1) | (index & 1);
                index >>= 1;
            }

            indices[i] = reversed;
        }
        return indices;
    }
    
    void InitializeButterflyTexture(int N, int[] bitReversedIndices)
    {
        int stages = Mathf.RoundToInt(Mathf.Log(N, 2)); // Total Stage (log_2(N))
        butterflyTexture = new Texture2D(stages, N, TextureFormat.RGBAFloat, false);
        
        for (int y = 0; y < butterflyTexture.height; y++)
        {
            for (int x = 0; x < butterflyTexture.width; x++)
            {
                butterflyTexture.SetPixel(x, y, new Color(0, 0, 0, 0)); // 初始化为全黑
            }
        }
        butterflyTexture.Apply();

        for (int x = 0; x < stages; x++)
        {
            int butterflySpan = 1 << x; // current span

            for (int y = 0; y < N; y++)
            {
                // twiddle factor
                // Equation (4.6)
                float k = (y * (float)N / (1 << (x + 1))) % N;
                float twiddleReal = Mathf.Cos(2.0f * Mathf.PI * k / N);
                float twiddleImag = Mathf.Sin(2.0f * Mathf.PI * k / N);

                bool butterflyWing = ((y % (1 << (x + 1))) < butterflySpan); // wing mark

                // int flippedStage = stages - x - 1;
                // first stage, bit reversed indices
                if (x == 0)
                {
                    if (butterflyWing)
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            bitReversedIndices[y],
                            bitReversedIndices[y + 1]
                        ));
                    }
                    else
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            bitReversedIndices[y - 1],
                            bitReversedIndices[y]
                        ));
                    }
                }
                else
                {
                    if (butterflyWing)
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            y,
                            y + butterflySpan
                        ));
                    }
                    else
                    {
                        butterflyTexture.SetPixel(x, y, new Color(
                            twiddleReal,
                            twiddleImag,
                            y - butterflySpan,
                            y
                        ));
                    }
                }
            }
        }
        butterflyTexture.Apply();
    }
    
    private Texture2D GenerateTexture(Complex[,] data)
    {
        Texture2D texture = new Texture2D(N, N, TextureFormat.RGFloat, false);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                float real = (float)data[i, j].Real;
                float imag = (float)data[i, j].Imaginary;
                // The real part of h0(k) is stored in the red channel
                // while the imaginary part is stored in the green channel
                texture.SetPixel(i, j, new Color(real, imag, 0));
            }
        }

        texture.Apply();
        return texture;
    }
    
    #endregion
    
    #region Phase 2 Height Amplitude Fourier Components
    Texture2D UpdateFourierComponents(Complex[,] h0k, Complex[,] h0MinusK, float time)
    {
        Texture2D hktTexture = new Texture2D(N, N, TextureFormat.RGFloat, false);

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                // Wave Vector \bold{k}
                // Equation (4.2)
                float kx = (2 * Mathf.PI * i - Mathf.PI * N ) / L;
                float kz = (2 * Mathf.PI * j - Mathf.PI * N ) / L;
                Vector2 k = new Vector2(kx, kz);

                float kLength = k.magnitude;
                if (kLength < 1e-6f) kLength = 1e-6f;

                // frequency omega(k)
                // Equation (3.7)
                float omega = Mathf.Sqrt(9.8f * kLength);

                // Fourier Component
                // Equation (3.8)
                Complex expPositive = h0k[i, j] * Complex.Exp(new Complex(0, omega * time));
                Complex expNegative = h0MinusK[i, j] * Complex.Exp(new Complex(0, -omega * time));
                Complex hkt = expPositive + expNegative;

                // Write into Texture
                hktTexture.SetPixel(i, j, new Color((float)hkt.Real, (float)hkt.Imaginary, 0));
            }
        }

        hktTexture.Apply();
        return hktTexture;
    }
    
    #endregion
    
    #region Phase 3 Horitontal 1D FFT
    Complex[,] ExtractDataFromTexture(Texture2D texture)
    {
        int N = texture.width;
        Complex[,] data = new Complex[N, N];

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                Color color = texture.GetPixel(i, j);
                data[i, j] = new Complex(color.r, color.g); 
            }
        }

        return data;
    }
    
   
    
    #endregion
    
    #region Phase 4 Vertical 1D FFT
    
    #endregion
    
    #region Phase 5 Inversion and Permutation
    
    #endregion
} 

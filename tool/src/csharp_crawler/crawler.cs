using System.Diagnostics;
using System.Globalization;
using System.IO.Compression;
using System.Net;
using System.Net.Http;
using System.Net.WebSockets;
using System.Text;
using System.Text.Json;

// ============================================================================
//  ICrawler —— 对外契约
// ============================================================================
public interface ICrawler
{
    /// <summary>取 Config 索引并完成无头浏览器初始化。</summary>
    void Init();

    /// <summary>按 ID 拼 URL 下载到 SAVEPATH。</summary>
    void GetOne(string ID);
}

// ============================================================================
//  Hasher —— 键字符串 -> UInt64（FNV-1a 64）
// ============================================================================
public static class Hasher
{
    private const ulong Offset = 14695981039346656037UL;
    private const ulong Prime = 1099511628211UL;

    public static ulong Of(string s)
    {
        if (string.IsNullOrEmpty(s)) return 0;
        ulong h = Offset;
        foreach (byte b in Encoding.UTF8.GetBytes(s))
        {
            h ^= b;
            h *= Prime;
        }
        return h;
    }
}

// ============================================================================
//  Stats —— 跨线程统计 + 控制台进度
// ============================================================================
public sealed class Stats
{
    private long _ok;
    private long _fail;
    private long _bytes;
    private long _done;
    private readonly object _console = new();
    private int _lineLen;
    private string _lastLive = "";

    public int Total { get; set; }
    public long Ok => Interlocked.Read(ref _ok);
    public long Fail => Interlocked.Read(ref _fail);
    public long Bytes => Interlocked.Read(ref _bytes);

    private void ClearLineLocked()
    {
        if (_lineLen > 0 && !Console.IsOutputRedirected)
            Console.Write("\r" + new string(' ', _lineLen) + "\r");
        _lineLen = 0;
        _lastLive = "";
    }

    private void Emit(string line)
    {
        lock (_console)
        {
            ClearLineLocked();
            Console.WriteLine(line);
        }
    }

    public void AddOk(string id, long bytes, long ms)
    {
        Interlocked.Increment(ref _ok);
        Interlocked.Add(ref _bytes, bytes);
        int done = (int)Interlocked.Increment(ref _done);
        Emit($"[{done,4}/{Total}] OK   {id,-5} {bytes / 1024.0,8:F1} KB {ms,6} ms");
    }

    public void AddFail(string id, string reason)
    {
        Interlocked.Increment(ref _fail);
        int done = (int)Interlocked.Increment(ref _done);
        Emit($"[{done,4}/{Total}] FAIL {id,-5} {reason}");
    }

    public void AddSkip(string id, string reason)
    {
        Interlocked.Increment(ref _ok);
        int done = (int)Interlocked.Increment(ref _done);
        Emit($"[{done,4}/{Total}] SKIP {id,-5} {reason}");
    }

    public void PrintLive()
    {
        if (Console.IsOutputRedirected) return;
        lock (_console)
        {
            int done = (int)(Ok + Fail);
            int pct = Total <= 0 ? 0 : (int)(done * 100L / Total);
            string s = $"下载中 {done}/{Total} ({pct,3}%)  ok={Ok} fail={Fail}  {Bytes / 1048576.0,7:F2} MB";
            if (s == _lastLive) return;
            Console.Write("\r" + s);
            _lineLen = s.Length;
            _lastLive = s;
        }
    }

    public void EndLive()
    {
        lock (_console) ClearLineLocked();
    }
}

// ============================================================================
//  Browser —— 无头 Chromium（CDP over WebSocket），全进程共享一个实例
//  仅用于初始化探测与 HTTP 失败时的兜底，不影响主下载路径。
// ============================================================================
public static class Browser
{
    private static readonly object Gate = new();
    private static Process _proc;
    private static ClientWebSocket _ws;
    private static int _nextId;
    private static bool _tried;

    public static bool Ready { get; private set; }
    public static string Version { get; private set; } = "-";

    public static void EnsureStarted(string chromePath, bool enable, int timeoutMs)
    {
        if (!enable)
        {
            Console.WriteLine("[browser] USE_BROWSER=false，跳过无头浏览器初始化");
            return;
        }
        lock (Gate)
        {
            if (_tried) return;
            _tried = true;
            try
            {
                Start(chromePath, timeoutMs);
                Ready = true;
                Console.WriteLine($"[browser] 无头浏览器就绪: {Version}");
            }
            catch (Exception ex)
            {
                Ready = false;
                Console.WriteLine($"[browser] 初始化失败，退回纯 HTTP: {ex.GetType().Name}: {ex.Message}");
            }
        }
    }

    /// <summary>CHROME_PATH 没配就按常见位置找 Chrome / Edge。</summary>
    private static string ResolveChrome(string configured)
    {
        if (!string.IsNullOrWhiteSpace(configured) && File.Exists(configured)) return configured;
        var candidates = new[]
        {
            Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
                "Google", "Chrome", "Application", "chrome.exe"),
            @"C:\Program Files\Google\Chrome\Application\chrome.exe",
            @"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
            @"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        };
        foreach (var c in candidates)
            if (File.Exists(c)) return c;
        return null;
    }

    private static void Start(string chromePath, int timeoutMs)
    {
        chromePath = ResolveChrome(chromePath);
        if (chromePath == null)
            throw new FileNotFoundException("找不到 Chromium/Chrome，可在 config.json 里写 CHROME_PATH");

        string userData = Path.Combine(Path.GetTempPath(), "crawler-chrome-" + Environment.ProcessId);
        if (Directory.Exists(userData)) Directory.Delete(userData, true);
        Directory.CreateDirectory(userData);

        var psi = new ProcessStartInfo(chromePath)
        {
            UseShellExecute = false,
            CreateNoWindow = true,
            // 必须重定向：否则 Chrome 会继承本进程 stdout 句柄，
            // 它活多久父进程的输出管道就关不掉（后台运行会表现为「卡住不退出」）
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            RedirectStandardInput = true,
        };
        foreach (var a in new[]
                 {
                     "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
                     "--disable-extensions", "--mute-audio", "--disable-background-networking",
                     "--disable-breakpad", "--no-crashpad",
                     "--remote-debugging-port=0", "--user-data-dir=" + userData, "about:blank"
                 })
            psi.ArgumentList.Add(a);

        _proc = Process.Start(psi) ?? throw new InvalidOperationException("chrome 启动失败");
        _proc.OutputDataReceived += (_, _) => { };
        _proc.ErrorDataReceived += (_, _) => { };
        _proc.BeginOutputReadLine();
        _proc.BeginErrorReadLine();
        AppDomain.CurrentDomain.ProcessExit += OnProcessExit;

        string portFile = Path.Combine(userData, "DevToolsActivePort");
        var sw = Stopwatch.StartNew();
        while (!File.Exists(portFile))
        {
            if (sw.ElapsedMilliseconds > timeoutMs) throw new TimeoutException("等待 DevToolsActivePort 超时");
            if (_proc.HasExited) throw new InvalidOperationException("chrome 提前退出，exit=" + _proc.ExitCode);
            Thread.Sleep(100);
        }

        string[] lines = Array.Empty<string>();
        for (int i = 0; i < 20; i++)
        {
            try
            {
                lines = File.ReadAllLines(portFile);
                if (lines.Length > 0 && lines[0].Trim().Length > 0) break;
            }
            catch (IOException) { }
            Thread.Sleep(100);
        }
        int port = int.Parse(lines[0].Trim(), CultureInfo.InvariantCulture);

        string wsUrl = null;
        using (var http = new HttpClient { Timeout = TimeSpan.FromSeconds(10) })
        {
            string json = http.GetStringAsync($"http://127.0.0.1:{port}/json/list").GetAwaiter().GetResult();
            using var doc = JsonDocument.Parse(json);
            foreach (var t in doc.RootElement.EnumerateArray())
            {
                if (t.TryGetProperty("type", out var type) && type.GetString() == "page" &&
                    t.TryGetProperty("webSocketDebuggerUrl", out var ws))
                {
                    wsUrl = ws.GetString();
                    break;
                }
            }
        }
        if (wsUrl == null) throw new InvalidOperationException("找不到 page target");

        _ws = new ClientWebSocket();
        _ws.Options.KeepAliveInterval = TimeSpan.FromSeconds(20);
        _ws.ConnectAsync(new Uri(wsUrl), CancellationToken.None).GetAwaiter().GetResult();

        Version = Evaluate("navigator.userAgent", 10000) ?? "unknown";
    }

    /// <summary>在页面上下文里执行 JS 表达式，返回字符串结果。</summary>
    public static string Evaluate(string expression, int timeoutMs)
    {
        lock (Gate)
        {
            if (_ws == null || _ws.State != WebSocketState.Open) throw new InvalidOperationException("CDP 未连接");
            int id = ++_nextId;
            string payload = "{\"id\":" + id +
                             ",\"method\":\"Runtime.evaluate\",\"params\":{\"expression\":" +
                             JsonSerializer.Serialize(expression) +
                             ",\"returnByValue\":true,\"awaitPromise\":true}}";
            var bytes = Encoding.UTF8.GetBytes(payload);
            _ws.SendAsync(new ArraySegment<byte>(bytes), WebSocketMessageType.Text, true, CancellationToken.None)
                .GetAwaiter().GetResult();

            var buf = new byte[1 << 16];
            var deadline = DateTime.UtcNow.AddMilliseconds(timeoutMs);
            while (DateTime.UtcNow < deadline)
            {
                var sb = new StringBuilder();
                WebSocketReceiveResult r;
                do
                {
                    r = _ws.ReceiveAsync(new ArraySegment<byte>(buf), CancellationToken.None).GetAwaiter().GetResult();
                    sb.Append(Encoding.UTF8.GetString(buf, 0, r.Count));
                } while (!r.EndOfMessage);

                using var doc = JsonDocument.Parse(sb.ToString());
                if (doc.RootElement.TryGetProperty("id", out var idEl) && idEl.GetInt32() == id)
                {
                    var root = doc.RootElement;
                    if (root.TryGetProperty("result", out var res) &&
                        res.TryGetProperty("result", out var inner) &&
                        inner.TryGetProperty("value", out var v))
                        return v.ValueKind == JsonValueKind.String ? v.GetString() : v.ToString();
                    return null;
                }
            }
            throw new TimeoutException("CDP evaluate 超时");
        }
    }

    /// <summary>把页面导航到某地址并读回纯文本（.pdb 是 text/plain，兜底可用）。</summary>
    public static string FetchText(string url, int timeoutMs)
    {
        lock (Gate)
        {
            int id = ++_nextId;
            string payload = "{\"id\":" + id + ",\"method\":\"Page.navigate\",\"params\":{\"url\":" +
                             JsonSerializer.Serialize(url) + "}}";
            var bytes = Encoding.UTF8.GetBytes(payload);
            _ws.SendAsync(new ArraySegment<byte>(bytes), WebSocketMessageType.Text, true, CancellationToken.None)
                .GetAwaiter().GetResult();
        }
        Thread.Sleep(1500);
        return Evaluate("document.body ? document.body.innerText : ''", timeoutMs);
    }

    /// <summary>进程退出（含 Ctrl+C）时的兜底清理，不加锁以免和主线程互等。</summary>
    private static void OnProcessExit(object sender, EventArgs e)
    {
        try { if (_proc is { HasExited: false }) _proc.Kill(true); } catch { }
    }

    public static void Dispose()
    {
        lock (Gate)
        {
            try { _ws?.CloseAsync(WebSocketCloseStatus.NormalClosure, "bye", CancellationToken.None).Wait(2000); } catch { }
            try { _ws?.Dispose(); } catch { }
            _ws = null;
            try
            {
                if (_proc is { HasExited: false })
                {
                    _proc.Kill(true);          // 连带子进程一起收掉
                    _proc.WaitForExit(5000);
                }
            }
            catch { }
            try { _proc?.Dispose(); } catch { }
            _proc = null;
            Ready = false;
        }
    }
}

// ============================================================================
//  Crawler —— 一个线程一个实例
// ============================================================================
public sealed class Crawler : ICrawler, IDisposable
{
    private readonly Config _cfg;
    private readonly Stats _stats;
    private readonly int _no;

    // Config 索引（键字符串的 UInt64 哈希码）
    private readonly ulong _kHead = Hasher.Of("API_HEAD");
    private readonly ulong _kEnd = Hasher.Of("API_END");
    private readonly ulong _kSave = Hasher.Of("SAVEPATH");
    private readonly ulong _kChrome = Hasher.Of("CHROME_PATH");
    private readonly ulong _kUseBrowser = Hasher.Of("USE_BROWSER");
    private readonly ulong _kTimeout = Hasher.Of("TIMEOUT_SEC");
    private readonly ulong _kRetry = Hasher.Of("RETRY");

    private string _apiHead;
    private string _apiEnd;
    private string _savePath;
    private string _chromePath;
    private bool _useBrowser;
    private int _timeoutSec;
    private int _retry;
    private bool _gz;

    private HttpClient _http;
    private bool _inited;

    public Crawler(Config cfg, int no, Stats stats)
    {
        _cfg = cfg;
        _no = no;
        _stats = stats;
    }

    public void Init()
    {
        // --- 1. 获得 Config 索引 ---
        _apiHead = _cfg.Get(_kHead, "https://files.rcsb.org/download/");
        _apiEnd = _cfg.Get(_kEnd, ".pdb");
        _savePath = _cfg.Get(_kSave, "./data/pdb");
        _chromePath = _cfg.Get(_kChrome, "");
        _useBrowser = _cfg.GetBool(_kUseBrowser, true);
        _timeoutSec = _cfg.GetInt(_kTimeout, 60);
        _retry = Math.Max(1, _cfg.GetInt(_kRetry, 3));

        // 容错：API_HEAD 允许不带 scheme、允许不写结尾的 /
        if (!_apiHead.StartsWith("http", StringComparison.OrdinalIgnoreCase)) _apiHead = "https://" + _apiHead;
        if (!_apiHead.EndsWith("/", StringComparison.Ordinal)) _apiHead += "/";

        _gz = _apiEnd.EndsWith(".gz", StringComparison.OrdinalIgnoreCase);

        Directory.CreateDirectory(_savePath);

        // --- 2. 完成无头浏览器初始化（全进程共享，只有第一个 Crawler 真正拉起）---
        Browser.EnsureStarted(_chromePath, _useBrowser, _timeoutSec * 1000);

        // --- 3. HTTP 会话 ---
        _http = new HttpClient(new SocketsHttpHandler
        {
            MaxConnectionsPerServer = 16,
            PooledConnectionLifetime = TimeSpan.FromMinutes(5),
        })
        {
            Timeout = TimeSpan.FromSeconds(_timeoutSec),
        };
        _http.DefaultRequestHeaders.TryAddWithoutValidation(
            "User-Agent",
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) " +
            "Chrome/153.0.0.0 Safari/537.36 rcsb-pdb-crawler/1.0");
        _http.DefaultRequestHeaders.TryAddWithoutValidation("Accept", "*/*");

        _inited = true;
    }

    public void GetOne(string ID)
    {
        if (!_inited) Init();

        string url = _apiHead + ID + _apiEnd;          // 拼 API_HEAD + ID + API_END
        string final = Path.Combine(_savePath, ID + ".pdb");
        string part = final + ".part";

        if (File.Exists(final) && LooksLikePdb(final))
        {
            _stats.AddSkip(ID, "已存在");
            return;
        }

        var sw = Stopwatch.StartNew();
        string lastError = null;

        for (int attempt = 1; attempt <= _retry; attempt++)
        {
            try
            {
                DownloadTo(url, part);
                long size = new FileInfo(part).Length;
                if (!LooksLikePdb(part)) throw new InvalidDataException($"内容不像 PDB（{size} B）");
                File.Move(part, final, true);
                _stats.AddOk(ID, size, sw.ElapsedMilliseconds);
                return;
            }
            catch (Exception ex)
            {
                lastError = ex.GetType().Name + ": " + ex.Message;
                TryDelete(part);
                if (attempt < _retry) Thread.Sleep(300 * attempt);
            }
        }

        // HTTP 连续失败 -> 用无头浏览器兜底（仅未压缩的纯文本情形）
        if (Browser.Ready && !_gz)
        {
            try
            {
                string text = Browser.FetchText(url, _timeoutSec * 1000);
                if (!string.IsNullOrWhiteSpace(text) && text.Contains("END"))
                {
                    File.WriteAllText(part, text, new UTF8Encoding(false));
                    File.Move(part, final, true);
                    _stats.AddOk(ID, new FileInfo(final).Length, sw.ElapsedMilliseconds);
                    return;
                }
            }
            catch (Exception ex)
            {
                lastError += " | browser: " + ex.Message;
            }
        }

        _stats.AddFail(ID, lastError ?? "unknown");
    }

    private void DownloadTo(string url, string path)
    {
        using var resp = _http.GetAsync(url, HttpCompletionOption.ResponseHeadersRead).GetAwaiter().GetResult();
        if (!resp.IsSuccessStatusCode)
            throw new HttpRequestException($"HTTP {(int)resp.StatusCode} {resp.ReasonPhrase}");

        using var src = resp.Content.ReadAsStream();
        using var dst = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 64 * 1024);
        if (_gz)
        {
            using var gz = new GZipStream(src, CompressionMode.Decompress);
            gz.CopyTo(dst, 64 * 1024);
        }
        else
        {
            src.CopyTo(dst, 64 * 1024);
        }
    }

    private static bool LooksLikePdb(string path)
    {
        try
        {
            var fi = new FileInfo(path);
            if (!fi.Exists || fi.Length < 200) return false;

            using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            var head = new byte[6];
            if (fs.Read(head, 0, 6) != 6) return false;
            string magic = Encoding.ASCII.GetString(head);
            string[] ok = { "HEADER", "ATOM  ", "HETATM", "TITLE ", "REMARK", "CRYST1", "EXPDTA", "MODEL ", "SEQRES" };
            bool good = false;
            foreach (var m in ok)
                if (magic == m) { good = true; break; }
            if (!good) return false;

            long tailLen = Math.Min(256, fi.Length);
            fs.Seek(-tailLen, SeekOrigin.End);
            var tail = new byte[tailLen];
            int n = fs.Read(tail, 0, (int)tailLen);
            return Encoding.ASCII.GetString(tail, 0, n).Contains("END");
        }
        catch
        {
            return false;
        }
    }

    private static void TryDelete(string path)
    {
        try { if (File.Exists(path)) File.Delete(path); } catch { }
    }

    public void Dispose()
    {
        try { _http?.Dispose(); } catch { }
        _http = null;
        _inited = false;
    }

    public override string ToString() => $"Crawler#{_no}";
}

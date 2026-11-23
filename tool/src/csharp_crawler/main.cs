using System.Diagnostics;
using System.Text;
using System.Text.Json;

// ============================================================================
//  Config —— 键的 UInt64 哈希码 -> 原始字符串值
// ============================================================================
public struct Config
{
    public Dictionary<ulong, string> Map;

    public int Count => Map?.Count ?? 0;

    public IEnumerable<KeyValuePair<ulong, string>> Items =>
        Map ?? Enumerable.Empty<KeyValuePair<ulong, string>>();

    public string Get(ulong key)
    {
        if (Map != null && Map.TryGetValue(key, out var v)) return v;
        return null;
    }

    public string Get(ulong key, string fallback) => Get(key) ?? fallback;

    public int GetInt(ulong key, int fallback) =>
        int.TryParse(Get(key), out var v) ? v : fallback;

    public bool GetBool(ulong key, bool fallback)
    {
        string s = Get(key);
        if (s == null) return fallback;
        s = s.Trim().ToLowerInvariant();
        return s is "1" or "true" or "yes" or "on";
    }
}

// ============================================================================
//  Main —— 全流程：读配置 -> 建 Crawler -> 抽 ID -> 并发下载 -> 收尾
// ============================================================================
public class Main
{
    private readonly ulong _kCrawlerNum = Hasher.Of("CRAWLER_NUM");
    private readonly ulong _kIdPath = Hasher.Of("ID_PATH");
    private readonly ulong _kIdNum = Hasher.Of("ID_NUM");

    public int ExitCode { get; private set; }

    public Main(string[] args)
    {
        Console.OutputEncoding = new UTF8Encoding(false);
        Console.InputEncoding = Encoding.UTF8;

        var crawlers = new List<Crawler>();
        string[] ids = Array.Empty<string>();

        try
        {
            // ---- 1. 解析当前目录 config.json ----
            string configPath = args.Length > 0
                ? args[0]
                : Path.Combine(Environment.CurrentDirectory, "config.json");

            Console.WriteLine("== RCSB PDB Crawler ==");
            Console.WriteLine($"工作目录 : {Environment.CurrentDirectory}");
            Console.WriteLine($"配置文件 : {configPath}");
            if (!File.Exists(configPath))
            {
                Console.WriteLine($"[错误] 找不到配置文件: {configPath}");
                ExitCode = 2;
                return;
            }

            // ---- 2. new Config，键哈希化成 UInt64 后统一用它索引 ----
            var cfg = new Config { Map = new Dictionary<ulong, string>() };
            var display = new List<(string Name, string Value)>();
            string raw = File.ReadAllText(configPath, Encoding.UTF8);
            using (var doc = JsonDocument.Parse(NormalizeBackslashes(raw), new JsonDocumentOptions
                   {
                       AllowTrailingCommas = true,
                       CommentHandling = JsonCommentHandling.Skip,
                   }))
            {
                if (doc.RootElement.ValueKind != JsonValueKind.Object)
                    throw new InvalidDataException("config.json 顶层必须是对象");
                Flatten(doc.RootElement, "", cfg.Map, display);
            }

            Console.WriteLine($"配置项   : {display.Count}");
            foreach (var d in display)
                Console.WriteLine($"  {Hasher.Of(d.Name),20}  {d.Name,-22} -> {d.Value}");

            int crawlerNum = Math.Max(1, cfg.GetInt(_kCrawlerNum, 4));
            string idPath = cfg.Get(_kIdPath, "./data/ids/entries.txt");
            int idNum = Math.Max(1, cfg.GetInt(_kIdNum, 20));

            // ---- 3. 按 CRAWLER_NUM new 出多个 Crawler 并初始化 ----
            var stats = new Stats();
            Console.WriteLine($"\n[1/5] 创建 {crawlerNum} 个 Crawler …");
            for (int i = 0; i < crawlerNum; i++)
            {
                var c = new Crawler(cfg, i + 1, stats);
                c.Init();                 // 取 Config 索引 + 无头浏览器初始化
                crawlers.Add(c);
            }
            Console.WriteLine($"[1/5] 完成，{crawlers.Count} 个 Crawler 就绪");

            // ---- 4. 从 ID_PATH 随机取 ID_NUM 个 ----
            Console.WriteLine($"\n[2/5] 读取 ID 池: {idPath}");
            if (!File.Exists(idPath))
                throw new FileNotFoundException("ID_PATH 指向的文件不存在", idPath);

            var pool = LoadPool(idPath);
            if (pool.Length == 0) throw new InvalidDataException("ID 池为空");

            int pick = Math.Min(idNum, pool.Length);
            ids = Sample(pool, pick);
            Console.WriteLine($"[2/5] 池内 {pool.Length} 条，随机抽取 {pick} 条（种子已固定，可复现）");

            // ---- 5. 并发下载 ----
            stats.Total = pick;
            Console.WriteLine($"\n[3/5] 开始下载 …\n");

            using var stop = new CancellationTokenSource();
            var progress = new Thread(() =>
            {
                while (!stop.IsCancellationRequested)
                {
                    stats.PrintLive();
                    Thread.Sleep(250);
                }
            }) { IsBackground = true };
            progress.Start();

            int cursor = -1;
            var sw = Stopwatch.StartNew();
            var workers = crawlers.Select(c => Task.Run(() =>
            {
                while (true)
                {
                    int i = Interlocked.Increment(ref cursor);
                    if (i >= ids.Length) break;
                    try { c.GetOne(ids[i]); }
                    catch (Exception ex) { stats.AddFail(ids[i], "异常: " + ex.Message); }
                }
            })).ToArray();

            Task.WaitAll(workers);
            sw.Stop();
            stop.Cancel();
            progress.Join();
            stats.EndLive();

            Console.WriteLine($"\n[4/5] 下载结束: 成功 {stats.Ok} / 失败 {stats.Fail} / 共 {pick}");
            Console.WriteLine($"      耗时 {sw.Elapsed.TotalSeconds:F1} s，落盘 {stats.Bytes / 1048576.0:F2} MB");
            ExitCode = stats.Fail == 0 ? 0 : 1;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"\n[致命] {ex.GetType().Name}: {ex.Message}");
            ExitCode = 3;
        }
        finally
        {
            // ---- 6. 释放内存 ----
            Console.WriteLine("[5/5] 释放资源 …");
            foreach (var c in crawlers)
            {
                try { c.Dispose(); } catch { }
            }
            crawlers.Clear();
            Browser.Dispose();
            ids = Array.Empty<string>();
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
            Console.WriteLine("[5/5] 完成");
        }
    }

    /// <summary>
    /// 把 JSON 递归拍平：嵌套对象用 "父.子" 作为键名，同时给叶子名留一个别名，
    /// 这样 {"PBD_API":{"API_HEAD":...}} 和 {"API_HEAD":...} 两种写法都能用同一套哈希索引。
    /// </summary>
    private static void Flatten(JsonElement obj, string prefix, Dictionary<ulong, string> map,
        List<(string Name, string Value)> display)
    {
        foreach (var prop in obj.EnumerateObject())
        {
            string path = string.IsNullOrEmpty(prefix) ? prop.Name : prefix + "." + prop.Name;
            if (prop.Value.ValueKind == JsonValueKind.Object)
            {
                Flatten(prop.Value, path, map, display);
                continue;
            }
            string value = prop.Value.ValueKind == JsonValueKind.String
                ? prop.Value.GetString()
                : prop.Value.GetRawText();
            map[Hasher.Of(path)] = value;
            display.Add((path, value));
            if (!string.IsNullOrEmpty(prefix)) map.TryAdd(Hasher.Of(prop.Name), value);
        }
    }

    /// <summary>
    /// 手写配置里的 Windows 路径通常只写一个反斜杠（"D:\a\b"），严格 JSON 会解析失败
    /// 或被 \r \t 之类悄悄吃掉。这里只把 \" \\ \/ 当转义，其余反斜杠一律当路径分隔符。
    /// </summary>
    private static string NormalizeBackslashes(string raw)
    {
        var sb = new StringBuilder(raw.Length + 32);
        bool inString = false;
        for (int i = 0; i < raw.Length; i++)
        {
            char c = raw[i];
            if (!inString)
            {
                if (c == '"') inString = true;
                sb.Append(c);
                continue;
            }
            if (c == '"')
            {
                inString = false;
                sb.Append(c);
                continue;
            }
            if (c == '\\')
            {
                char n = i + 1 < raw.Length ? raw[i + 1] : '\0';
                if (n is '"' or '\\' or '/')
                {
                    sb.Append(c).Append(n);
                    i++;
                }
                else
                {
                    sb.Append("\\\\");
                }
                continue;
            }
            sb.Append(c);
        }
        return sb.ToString();
    }

    /// <summary>ID 池：既支持 JSON 数组（RCSB holdings 原样），也支持一行一个。</summary>
    private static string[] LoadPool(string path)
    {
        string raw = File.ReadAllText(path, Encoding.UTF8).Trim();
        if (raw.StartsWith("["))
        {
            using var doc = JsonDocument.Parse(raw);
            return doc.RootElement.EnumerateArray()
                .Select(e => e.GetString())
                .Where(s => !string.IsNullOrWhiteSpace(s))
                .Select(s => s.Trim())
                .ToArray();
        }
        return raw.Split('\n')
            .Select(l => l.Trim())
            .Where(l => l.Length > 0)
            .ToArray();
    }

    /// <summary>部分 Fisher–Yates，固定种子便于复现。</summary>
    private static string[] Sample(string[] pool, int count)
    {
        var rnd = new Random(20260101);
        var idx = new int[pool.Length];
        for (int i = 0; i < idx.Length; i++) idx[i] = i;
        for (int i = 0; i < count; i++)
        {
            int j = rnd.Next(i, idx.Length);
            (idx[i], idx[j]) = (idx[j], idx[i]);
        }
        var result = new string[count];
        for (int i = 0; i < count; i++) result[i] = pool[idx[i]];
        return result;
    }
}

public static class Program
{
    public static int Main(string[] args) => new Main(args).ExitCode;
}

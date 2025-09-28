# API Quota Monitor Tool

This tool helps you diagnose and troubleshoot API rate limiting issues with the Groq API used in Fastalysis.

## Usage

```bash
# Basic API status check
python api_quota_monitor.py --check

# Test rate limits with multiple requests
python api_quota_monitor.py --test --count 5

# Check environment variables
python api_quota_monitor.py --env

# Check service status
python api_quota_monitor.py --status
```

## Fixing Rate Limit Issues

If you're experiencing rate limit errors (`HTTP 429 Too Many Requests`), try these solutions:

1. **Add delays between requests** - Ensure your code waits at least 1 second between API calls
2. **Implement retry with exponential backoff** - Automatically retry failed requests with increasing delays
3. **Monitor usage patterns** - Use this tool to identify when rate limits are being hit
4. **Reduce concurrency** - Avoid multiple simultaneous requests to the API

## Environment Setup

Make sure your API keys are properly configured in `.env.local` file:

```
GROQ_API_KEY=gsk_your_key_here
```

## Interpreting Results

- **Success responses (HTTP 200)**: Your API key and quota are working correctly
- **Rate limit errors (HTTP 429)**: You've exceeded your allowed request rate
- **Authentication errors (HTTP 401)**: Your API key is invalid or missing
- **Other errors**: May indicate service issues or network problems

## Advanced Troubleshooting

For more comprehensive troubleshooting, see `RATE_LIMIT_GUIDE.md` in the project root.

For specific code fixes, see the `api_utils.py` module which implements rate limiting and retry logic.
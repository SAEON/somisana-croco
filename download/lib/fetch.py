import asyncio
import aiohttp

async def fetch(url):
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            return await response.text()

async def main():
    urls = [
        "http://example.com",
        "http://example.org",
        "http://example.net"
    ]
    tasks = []
    for url in urls:
        task = asyncio.ensure_future(fetch(url))
        tasks.append(task)
    responses = await asyncio.gather(*tasks)
    print(responses)

# loop = asyncio.get_event_loop()
# loop.run_until_complete(main())